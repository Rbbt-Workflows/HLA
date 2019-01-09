module HLA

  input :bam, :file, "Tumor BAM", nil, :nofile => true
  input :reference, :select, "Reference code", "hg19", :select_options => %w(b37 hg38)
  task :polysolver => :text do |bam,reference|
    directory = file('workdir').find
    bam = File.expand_path(bam)
    prepared_bam = Samtools.prepare_BAM(bam)
    bam_file = directory['file.bam']
    begin
      Open.ln_h File.expand_path(bam), bam_file
      Open.ln_h prepared_bam + '.bai', bam_file + '.bai'
      POLYSOLVER.run(directory, reference)
    ensure
      Open.rm bam_file if File.exists?(bam_file)
      Open.rm bam_file + '.bai' if File.exists?(bam_file + '.bai')
    end
    Open.read directory["winners.hla.txt"]
  end

  input :BAM, :file, "BAM file", nil, :nofile => true
  input :reference, :select, "Organism reference", 'b37', :select_options => %w(b37 hg38)
  task :SOAPHLA => :text do |bam, reference|
    dir = file('output')
    prepared_bam = Samtools.prepare_BAM(bam)
    SOAPHLA.run(prepared_bam, dir, reference)
    samples = dir.glob("*").collect{|f| File.basename(f)}
    files = samples.collect{|s| dir[s][s + '.type']}
    files.first.read
    #dir.type.read
  end

  input :files, :array, "List of fastq files", nil, :nofile => true
  dep :HLA_reads, :fastq1 => :placeholder, :fastq2 => nil do |jobname, options|
    options[:files].collect do |file|
      {:inputs => {:fastq1 => file}, :jobname => jobname}
    end
  end
  task :OptiType => :text do

    fastqs = dependencies.collect{|d| d.path }

    output = files_dir
    CMD.cmd("python #{Rbbt.software.opt.OptiType.produce["OptiTypePipeline.py"].find} -i #{fastqs.collect{|f| "'#{f}'"} * " "} --dna -v -o '#{output}' ")

    res = Dir.glob(files_dir + '/**/*.tsv').first

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
