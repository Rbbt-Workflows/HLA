module HLA

  input :fastq1, :file, "FASTQ file", nil, :nofile => true
  input :fastq2, :file, "FASTQ file 2", nil, :nofile => true
  input :filter_aligner, :select, "Aligner to use", :bwa, :select_options => %w(novoalign bwa razers3 bowtie)
  extension "fastq.gz"
  dep_task :HLA_reads, HTS, :bwa_filter, :reference => HLA_REFERENCE do |jobname,options|
    case options[:filter_aligner].to_s
    when "bwa"
      {:inputs => options, :workflow => HTS, :task => 'bwa_filter', :jobname => jobname}
    when "razers3"
      {:inputs => options, :workflow => HTS, :task => 'razers3_filter', :jobname => jobname}
    when "novoalign"
      {:inputs => options, :workflow => HTS, :task => 'novoalign_filter', :jobname => jobname}
    when "bowtie"
      {:inputs => options, :workflow => HTS, :task => 'bowtie_filter', :jobname => jobname}
    when ""
      raise ParameterException, "No aligner provided"
    else
      raise ParameterException, "Cannot find aligner: #{options[:aligner].to_s}"
    end
  end

end
