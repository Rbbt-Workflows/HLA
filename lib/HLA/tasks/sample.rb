Workflow.require_workflow "PVacSeq"

module Sample

  dep :BAM, :compute => [:produce, :canfail]
  dep :BAM_normal, :compute => [:produce, :canfail]
  dep HLA, :SOAPHLA, :compute => [:produce, :canfail], :BAM => :BAM_normal do |jobname,options,dependencies|
    options = add_sample_options jobname, options
    options[:BAM] = options[:bam] = :BAM if dependencies.flatten.reject{|dep| (dep.dependencies.empty? ) || (dep.error? && ! dep.recoverable_error?) }.select{|dep|  dep.task_name == :BAM_normal}.first.nil?
    {:inputs => options, :jobname => jobname}
  end
  dep HLA, :polysolver, :compute => [:produce, :canfail], :bam => :BAM_normal do |jobname,options,dependencies|
    options = add_sample_options jobname, options
    options[:BAM] = options[:bam] = :BAM if dependencies.flatten.reject{|dep| (dep.dependencies.empty? ) || (dep.error? && ! dep.recoverable_error?) }.select{|dep|  dep.task_name == :BAM_normal}.first.nil?
    {:inputs => options, :jobname => jobname.gsub(":", '_')}
  end
  dep HLA, :OptiType, :compute => [:produce, :canfail] do |sample,options|

    fastq_files = nil

    nsample = nil
    sample_study = Sample.sample_study(sample)
    sample_files = nil
    [sample, sample.sub('_normal', ''), sample + '_normal', [sample_study, "normal"] * ":"].uniq.each do |normal_sample|
      nsample = normal_sample
      sample_files = Sample.sample_files normal_sample if sample_study == Sample.sample_study(nsample)
      break if sample_files && fastq_files = sample_files[:FASTQ]
    end


    if sample_files && fastq_files = sample_files[:FASTQ]
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
  dep HLA, :arcasHLA, :compute => [:produce, :canfail], :bam => :BAM_normal do |jobname,options,dependencies|
    options = add_sample_options jobname, options
    options[:bam] = :BAM if dependencies.flatten.select{|dep| dep.task_name == :BAM_normal}.first.nil?
    {:inputs => options, :jobname => jobname.gsub(":", '_')}
  end
  task :hla_typing => :tsv do
    tsv = TSV.setup({}, :key_field => "Allele", :fields => ["Method"], :type => :flat)
    dependencies.each do |dep|
      next if dep.error?
      case dep.task_name
      when :BAM, :BAM_normal
        next
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
      else
        dep.path.read.split(/\s/).select{|w| w.include? "*"}.each do |allele|
          tsv[allele] ||= []
          tsv[allele] << dep.task_name
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
  dep_task :pvacseq, PVacSeq, :analysis, :vcf_file => :vcf_file, :alleles => :hla_genotype do |jobname,options,dependencies|
    options = add_sample_options jobname, options
    options.merge!(:organism => Organism.organism_for_build(options[:reference]))
    {:inputs => options, :jobname => jobname}
  end
  
  dep :pvacseq, :compute => :produce
  task :neo_epitopes => :tsv do

    organism = self.recursive_inputs[:organism] 

    parser = TSV::Parser.new step(:pvacseq).join, :type => :list
    tsv = TSV.setup({}, parser.options.merge(:key_field => "Genomic Mutation", :namespace => organism, :fields => parser.fields[4..-1]))
    TSV.traverse parser, :into => tsv do |chr, values|
      start, eend, ref, mut, *rest = values
      start = start.to_i
      start = start + 1 if ref.length == 1 && mut.length == 1
      pos, muts = Misc.correct_vcf_mutation start, ref, mut

      mutation = [chr, pos, muts.first] * ":"
      [mutation, values[4..-1]]
    end
    
    parser = TSV::Parser.new step(:pvacseq).file("output/MHC_Class_I/#{clean_name}.all_epitopes.tsv"), :header_hash => '', :type => :list
    all = TSV.setup({}, parser.options.merge(:key_field => "Genomic Mutation", :namespace => organism, :fields => parser.fields[4..-1]))
    TSV.traverse parser, :into => all do |chr, values|
      start, eend, ref, mut, *rest = values
      start = start.to_i
      start = start + 1 if ref.length == 1 && mut.length == 1
      pos, muts = Misc.correct_vcf_mutation start, ref, mut

      mutation = [chr, pos, muts.first] * ":"
      [mutation, values[4..-1]]
    end
    
    Open.write(file('all.tsv'), all.to_s)

    tsv
  end

  dep :mi
  dep_task :epitopes, HLA, :epitopes, :mutated_isoforms => :mi do |jobname, options|
    options = add_sample_options jobname, options
    options.merge!(:organism => Organism.organism_for_build(options[:reference])) if options[:organism].nil?
    {:inputs => options, :jobname => jobname}
  end

  dep :epitopes
  dep :hla_genotype
  dep_task :mhcFlurry, HLA, :mhcFlurry, "HLA#epitopes" => :epitopes, :alleles => :hla_genotype do |jobname, options|
    options = add_sample_options jobname, options
    options.merge!(:organism => Organism.organism_for_build(options[:reference])) if options[:organism].nil?
    {:inputs => options, :jobname => jobname}
  end

  dep :salmon
  task :gene_expression => :tsv do
    dumper = TSV::Dumper.new :key_field => "Ensembl Transcript ID", :fields => ["TPM"], :type => :single, :cast => :to_f
    dumper.init
    TSV.traverse step(:salmon), :header_hash => "", :type => :list, :into => dumper do |name,values|
      length, eff_length, tpm, reads = values
      transcript = name.split(".").first
      [transcript, tpm]
    end
  end

  dep :expanded_vcf
  dep :genomic_mutations
  dep :genomic_mutation_consequence
  dep :epitopes 
  dep :mhcFlurry 
  dep :gene_expression, :compute => :canfail
  dep :RNA_BAM, :compute => :canfail
  dep HTS, :BAM_position_pileup, :BAM => :RNA_BAM, :compute => :canfail, :positions => :genomic_mutations do |jobname,options,dependencies|
    if rna_bam = dependencies.select{|d| d.task_name == :RNA_BAM}.first
      reference = rna_bam.recursive_inputs[:reference]
      {:inputs => {:reference => reference}.merge(options), :jobname => jobname}
    else
      nil
    end
  end
  task :epitope_features => :tsv do
    organism = step(:genomic_mutation_consequence).recursive_inputs[:organism]
    organism = organism.load if Step === organism

    tsv = step(:epitopes).load.reorder "Mutated", nil, :zipped => true, :merge => true

    prot2enst = Organism.transcripts(organism).tsv :key_field => "Ensembl Protein ID", :fields => ["Ensembl Transcript ID"], :persist => true, :type => :single
    tsv = tsv.attach step(:genomic_mutation_consequence), :fields => ["Genomic Mutation"]
    tsv.add_field "Ensembl Transcript ID" do |mi,values|
      values["Mutated Isoform"].collect do |mi|
        protein = mi.split(":").first
        enst = prot2enst[protein]
        enst
      end
    end

    tsv = tsv.attach step(:gene_expression)
    tsv = tsv.attach step(:expanded_vcf)

    if step(:BAM_position_pileup).done?
      pileup = step(:BAM_position_pileup).load
      tsv.add_field "RNA Support" do |k,v|
        mutations = v["Genomic Mutation"]
        mutations.collect{|mutation|
          position = mutation.split(":").values_at(0,1) * ":"
          values = pileup[position]
          if values
            values["Alt count"]
          end
        }
      end
      tsv.add_field "RNA VAF" do |k,v|
        mutations = v["Genomic Mutation"]
        mutations.collect{|mutation|
          position = mutation.split(":").values_at(0,1) * ":"
          values = pileup[position]
          if values
            values["Alt count"].to_f / values["Ref count"].to_f
          end
        }
      end
    end

    mhcFlurry = step(:mhcFlurry).load
    mhcFlurry.add_field "Mutated" do |k,v|
      v["Mutated epitope"]
    end

    alleles = mhcFlurry.column("Allele").values.flatten.compact.uniq

    orig_fields = mhcFlurry.fields
    alleles.each do |allele|
      mhcFlurry.fields = orig_fields.collect{|f| f.include?("mhcflurry") ? f + " (#{allele})" : f}
      tsv = tsv.attach mhcFlurry.select("Allele" => allele), :fields => mhcFlurry.fields.select{|f| f.include? "mhcflurry"}
    end

    tsv
  end

end
