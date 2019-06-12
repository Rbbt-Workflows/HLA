Workflow.require_workflow "PVacSeq"

module Sample

  dep :BAM, :compute => [:produce, :canfail]
  dep :BAM_normal, :compute => [:produce, :canfail]
  dep HLA, :SOAPHLA, :compute => [:produce, :canfail], :BAM => :BAM_normal do |jobname,options,dependencies|
    options = add_sample_options jobname, options
    options[:BAM] = :BAM if dependencies.flatten.select{|dep| dep.done? && dep.task_name == :BAM_normal}.first.nil?
    {:inputs => options, :jobname => jobname}
  end
  dep HLA, :polysolver, :compute => [:produce, :canfail], :bam => :BAM_normal do |jobname,options,dependencies|
    options = add_sample_options jobname, options
    options[:bam] = :BAM if dependencies.flatten.select{|dep| dep.done? &&  dep.task_name == :BAM_normal}.first.nil?
    {:inputs => options, :jobname => jobname.gsub(":", '_')}
  end
  dep HLA, :OptiType, :compute => :canfail do |sample,options|

    fastq_files = nil

    nsample = nil
    sample_study = Sample.sample_study(sample)
    sample_files = nil
    [sample, sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
      break if sample_files && fastq_files = sample_files[:FASTQ]
    end

    if fastq_files = sample_files[:FASTQ]
      options = add_sample_options nsample, options
      options = options.merge({:files => fastq_files.flatten.uniq})
      {:inputs => options, :jobname => sample}
    else
      {}
    end
  end
  dep HLA, :xHLA, :compute => [:produce, :canfail], :bam => :BAM_normal do |jobname,options,dependencies|
    options = add_sample_options jobname, options
    options[:bam] = :BAM if dependencies.flatten.select{|dep| dep.task_name == :BAM_normal}.first.nil?
    {:inputs => options, :jobname => jobname.gsub(":", '_')}
  end
  task :hla_typing => :tsv do
    tsv = TSV.setup({}, :key_field => "Allele", :fields => ["Method"], :type => :flat)
    dependencies.each do |dep|
      next if dep.error?
      case dep.task_name
      when :xHLA
        dep.load.each do |allele|
          tsv[allele] ||= []
          tsv[allele] << "xHLA"
        end
      when :OptiType
        dep.path.read.split("\n").last.split("\t").select{|p| p.include? "*"}.each do |allele|
          tsv[allele] ||= []
          tsv[allele] << "OptiType"
        end
      when :SOAPHLA
        TSV.traverse dep, :type => :array do |line|
          allele1, allele2, *rest = line.split("\t")
          tsv[allele1] ||= []
          tsv[allele1] << "SOAPHLA"
        end
      when :polysolver
        TSV.traverse dep, :type => :array do |line|
          gene, allele1, allele2 = line.split("\t")
          tsv[allele1] ||= []
          tsv[allele2] ||= []
          tsv[allele1] << "Polysolver"
          tsv[allele2] << "Polysolver"
        end
      end
    end

    new = tsv.annotate({})
    tsv.each do |allele,dbs|
      normal = allele.upcase.sub("HLA",'').gsub("_", ":")
      normal.sub!(/^:/,'')
      normal.sub!(/([A-Z]):/,'\1*')
      next unless normal =~ /\d/
      new[normal] ||= []
      new[normal].concat dbs.uniq
    end
    new
  end

  dep :hla_typing do |sample,options,dependencies|
    nsample = nil
    sample_files = nil

    sample_study = Sample.sample_study(sample)
    [sample + '_normal', [sample_study, "normal"] * ":"].each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
      break if sample_files
    end

    if sample_files
      {:jobname => nsample, :inputs => options}
    else
      {:jobname => sample, :inputs => options}
    end
  end
  task :hla_genotype => :array do
    adbs = {}
    TSV.traverse step(:hla_typing) do |allele, dbs|
      allele = allele.first if Array ===  allele
      next if allele.split(":").length < 2
      short = allele.split(":")[0..1] * ":"
      adbs[short] ||= []
      adbs[short].concat dbs
    end
    adbs.keys
  end

  
  dep :vcf_file, :compute => :produce
  dep :hla_genotype, :compute => :produce
  dep PVacSeq, :analysis, :vcf_file => :vcf_file, :alleles => :hla_genotype do |jobname,options,dependencies|
    options = add_sample_options jobname, options
    options.merge!(:organism => Organism.organism_for_build(options[:reference]))
    {:inputs => options, :jobname => jobname}
  end
  task :neo_epitopes => :tsv do

    organism = self.recursive_inputs[:organism] 
    parser = TSV::Parser.new step(:analysis).join, :type => :list
    dumper = TSV::Dumper.new parser.options.merge(:key_field => "Genomic Mutation", :namespace => organism, :fields => parser.fields[4..-1])
    dumper.init
    TSV.traverse parser, :into => dumper do |chr, values|
      start, eend, ref, mut, *rest = values
      start = start.to_i
      start = start + 1 if ref.length == 1 && mut.length == 1
      pos, muts = Misc.correct_vcf_mutation start, ref, mut

      mutation = [chr, pos, muts.first] * ":"
      [mutation, values[4..-1]]
    end
  end


end
