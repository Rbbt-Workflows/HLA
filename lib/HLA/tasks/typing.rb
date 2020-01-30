module HLA

  input :bam, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38)
  task :polysolver => :text do |bam,reference|
    directory = file('workdir').find
    bam = File.expand_path(bam)
    prepared_bam = Samtools.prepare_BAM(bam)
    bam_file = directory['file.bam']
    begin
      Open.ln_h File.expand_path(bam), bam_file
      Open.ln_h prepared_bam + '.bai', bam_file + '.bai'
      POLYSOLVER.run(directory, reference)
    rescue Exception
      raise RbbtException, "Polysolver could not run"
    ensure
      Open.rm bam_file if File.exists?(bam_file)
      Open.rm bam_file + '.bai' if File.exists?(bam_file + '.bai')
    end
    Open.read directory["winners.hla.txt"]
  end

  input :bam, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "hg38", :select_options => %w(hg38)
  task :xHLA => :array do |bam,reference|
    raise ParameterException, "Only hg38 is allowed" if reference != "hg38"

    directory = file('workdir').find
    bam = File.expand_path(bam)
    prepared_bam = Samtools.prepare_BAM(bam)
    bam_file = directory['file.bam']
    begin
      Open.ln_h File.expand_path(bam), bam_file
      Open.ln_h prepared_bam + '.bai', bam_file + '.bai'
      XHLA.run(directory, reference)
    rescue Exception
      raise RbbtException, "XHLA could not run"
    ensure
      Open.rm bam_file if File.exists?(bam_file)
      Open.rm bam_file + '.bai' if File.exists?(bam_file + '.bai')
    end

    JSON.parse(file('workdir/report-sample-hla.json').read)["hla"]["alleles"]
  end

  input :BAM, :file, "BAM file", nil, :nofile => true
  input :reference, :select, "Organism reference", 'b37', :select_options => %w(b37)
  task :SOAPHLA => :text do |bam, reference|
    raise ParameterException, "Only b37 is allowed" if reference != "b37"
    dir = file('output')
    prepared_bam = Samtools.prepare_BAM(bam)
    SOAPHLA.run(prepared_bam, dir, reference)
    samples = dir.glob("*").collect{|f| File.basename(f)}
    files = samples.collect{|s| dir[s][s + '.type']}
    files.first.read
  end

  input :BAM, :file, "BAM file", nil, :nofile => true
  task :arcasHLA => :array do |bam|
    alleles = %w(A B C DPB1 DQB1 DQA1 DRB1)
    output = file('output')
    bam = Samtools.prepare_BAM bam
    ArcasHLA.run(bam,alleles,output).values.flatten.uniq
  end
  
  input :files, :array, "List of fastq files", nil, :nofile => true
  dep :HLA_reads, :fastq1 => :placeholder, :fastq2 => nil do |jobname, options|
    if options[:files] 
      options[:files].collect do |file|
        {:inputs => {:fastq1 => file}, :jobname => jobname}
      end
    else
      []
    end
  end
  task :arcasHLA_reads => :array do 
    alleles = %w(A B C DPB1 DQB1 DQA1 DRB1)
    output = file('output')
    inputs = file('input')

    fastqs = dependencies.collect{|dep| dep.path}
    fastqs_2 = dependencies.select{|d| d.recursive_inputs[:fastq1] =~ /_2\.fastq/}.collect{|d| d.path }
    fastqs_1 = fastqs - fastqs_2

    fastq1 = inputs["fastq1.fastq"]
    Open.open(fastq1, :mode => 'w') do |fin|
        fastqs_1.each do |fastq|
          Open.open(fastq) do |fout|
            Misc.consume_stream fout, false, fin, false
          end
        end
        fin.close
    end if fastqs_1.any?

    fastq2 = inputs["fastq2.fastq"]
    Open.open(fastq2, :mode => 'w') do |fin|
      fastqs_2.each do |fastq|
        Open.open(fastq) do |fout|
          Misc.consume_stream fout, false, fin, false
        end
      end
      fin.close
    end if fastqs_2.any?

    ArcasHLA.run_fasta(fastq1, fastq2,alleles,output).values.flatten.uniq
  end


  input :files, :array, "List of fastq files", nil, :nofile => true
  dep :HLA_reads, :fastq1 => :placeholder, :fastq2 => nil do |jobname, options|
    if options[:files] 
      files = options[:files]
      files = [files] unless Array === files
      files.collect do |file|
        file = file.path if Step === file
        {:inputs => {:fastq1 => file}, :jobname => jobname}
      end
    else
      []
    end
  end
  task :OptiType => :text do

    fastqs = dependencies.collect{|d| d.path }

    output = file('output')
    inputs = file('input')
    `chmod 777 -R #{output}`
    begin
      fastqs_2 = dependencies.select{|d| d.recursive_inputs[:fastq1] =~ /_2\.fastq/}.collect{|d| d.path }
      fastqs_1 = fastqs - fastqs_2

      fastq1 = inputs["fastq1.fastq"]
      Open.open(fastq1, :mode => 'w') do |fin|
        fastqs_1.each do |fastq|
          Open.open(fastq) do |fout|
            Misc.consume_stream fout, false, fin, false
          end
        end
        fin.close
      end if fastqs_1.any?

      fastq2 = inputs["fastq2.fastq"]
      Open.open(fastq2, :mode => 'w') do |fin|
        fastqs_2.each do |fastq|
          Open.open(fastq) do |fout|
            Misc.consume_stream fout, false, fin, false
          end
        end
        fin.close
      end if fastqs_2.any?

      fastqs = []
      fastqs << fastq1 if fastqs_1.any?
      fastqs << fastq2 if fastqs_2.any?

      #Docker.run('fred2/optitype',"-i #{fastqs.collect{|p| "/job/" << File.basename(p) } * " "} --dna -v -o /data/", :mounts => {"/data/" => output}, :directory => inputs)
      CMD.cmd_log("python #{Rbbt.software.opt.OptiType.produce["OptiTypePipeline.py"].find} -i #{fastqs.collect{|f| "'#{f}'"} * " "} --dna -v -o '#{output}'")
    rescue
      Log.exception $!
      raise RbbtException, "Error in python: #{$!.message}"
    ensure
      Open.rm_rf inputs
    end

    res = Dir.glob(output + '/**/*.tsv').first

    FileUtils.cp res, self.path
    nil
  end

  input :normalBAM, :file, "normal BAM file", nil, :nofile => true
  input :otherBAMs, :array, "other (tumor) BAM files", [], :nofile => true
  input :reference, :select, "Organism reference", 'b37', :select_options => %w(b37 hg38)
  input :alleles, :array, "HLA allele codes"
  input :purity, :float, "Tumor purity", 1.0
  input :ploidy, :float, "Tumor ploidy", 2.0
  task :LOHHLA => :text do |normalBAM, otherBAMs, reference, alleles,purity,ploidy|
    prepared_normal_bam = Samtools.prepare_BAM(normalBAM)
    prepared_other_bams = otherBAMs.collect{|bam| Samtools.prepare_BAM(bam) }

    LOHHLA.run(file('workdir'), alleles, prepared_normal_bam, prepared_other_bams, reference, purity, ploidy)
    "#" + Open.read(file('workdir').glob("*HLAlossPrediction_CI*").first)
  end
end
