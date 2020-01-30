module HLA

  helper :allele_info do |alleles|
    allele2code = IGMT.identifiers.index :target => "Code", :persist => true
    allele2size = IGMT.identifiers.tsv :fields => "Size", :type => :single, :persist => true

    allele_info = {}
    alleles.each do |allele|
      code = allele2code[allele]
      size = allele2size[allele]
      if code.nil?
        parts = allele.split(":")
        new_allele = parts.values_at(0,1,2,3).collect{|p|
          p || '01'
        } * ":"
        allele = allele + "_used_#{new_allele}"
        code = allele2code[new_allele]
        size = allele2size[new_allele]
      end
      next if code.nil?
      allele_info[allele] = [code, size]
    end

    allele_info
  end

  input :alleles, :array, "Alleles to use", nil
  extension 'fasta.gz'
  task :HLA_fasta => :binary do |alleles|
    alleles = Open.read(alleles).split("\n") if Misc.is_filename?(alleles)
    if alleles and alleles.any?
      allele_file = file('alleles.txt')
      allele_info = self.allele_info alleles
      codes = allele_info.collect{|k,v| v.first }.flatten.uniq
      Open.write(allele_file, codes*"\n")
      CMD.cmd_log("cat '#{HLA_REFERENCE.find}' | tr '\\n' '#' | sed 's/#>/\\n>/g;s/bp#/bp\\n/g'|grep -F -w -f '#{allele_file}' -A 1 | sed 's/#/\\n/g'| grep -v '^--$' | bgzip -c > #{self.tmp_path}")
    else
      CMD.cmd_log("bgzip #{HLA_REFERENCE} -c > '#{self.tmp_path}'")
    end
    nil
  end

  dep HLA, :HLA_fasta, :compute => :produce
  dep HLA, :HLA_reads, :reference => :HLA_fasta
  dep HTS, :BAM_bwa, :reference => :HLA_fasta, :fastq1 => :HLA_reads, :fastq2 => nil do |jobname,options|
    options[:fastq1] = :HLA_reads
    options[:fastq2] = nil
    case options[:aligner].to_s
    when "bwa"
      {:inputs => options.merge(:bwa_mem_args => '-a', :reference => :HLA_fasta), :task => :BAM_bwa, :workflow => HTS, :jobname => jobname}
    when "razers3"
      {:inputs => options.merge(:m => 1000, :reference => :HLA_fasta), :workflow => HTS, :task => 'razers3_BAM', :jobname => jobname}
    when "novoalign"
      {:inputs => options, :workflow => HTS, :task => 'novoalign_BAM', :jobname => jobname}
    when "bowtie"
      {:inputs => options.merge(:bowtie_args => '--end-to-end -a --no-unal'), :workflow => HTS, :task => 'bowtie_BAM', :jobname => jobname}
    when ""
      raise ParameterException, "No aligner provided"
    else
      raise ParameterException, "Cannot find aligner: #{options[:aligner].to_s}"
    end
  end
  input :aligner, :select, "Aligner to use", :bwa, :select_options => %w(bwa razers3 novoalign bowtie)
  extension :bam
  task :HLA_BAM => :binary do |aligner|
    if aligner.to_s != 'razers3' and aligner.to_s != 'novoalign'
      CMD.cmd_log("samtools view -h '#{dependencies.last.path}' | grep '^@\\|XM:i:0\\s' | samtools view -h -b   - > '#{self.tmp_path}'")
    else
      Open.ln_s dependencies.last.path, self.tmp_path
    end
    nil
  end


  dep :HLA_BAM, :compute => :produce
  input :alleles, :array, "HLA allele codes"
  input :image_format, :select, "Format for images", :png, :select_options => %w(png svg jpeg)
  input :complete_alleles, :boolean, "Add missing fields in allele names when not found", false
  task :allele_alignment => :array do |alleles,image_format,complete_alleles|

    hla_reference = step(:HLA_BAM).recursive_inputs[:reference]
    hla_reference = hla_reference.path if Step === hla_reference

    snapshot = file('snapshots')

    alignment = Samtools.prepare_BAM(step(:HLA_BAM).path)

    allele_info = self.allele_info alleles
    if Open.gzip?(hla_reference)
      Open.mkdir files_dir
      hla_reference_unzip = file('hla_reference.fasta')
      CMD.cmd("gunzip '#{hla_reference}' -c > '#{hla_reference_unzip}'")
      hla_reference = hla_reference_unzip
    end

    if allele_info.any?

      IGV.run <<-EOF, 10_000, 10_000
new
genome #{hla_reference}
load #{alignment}
snapshotDirectory #{snapshot}
maxPanelHeight 10000
#{
    allele_info.collect do |allele, info|
      code, size = info
      "goto #{code}:1-#{size}" << "\n" << "snapshot #{allele}.#{image_format}"
    end * "\n"
}
exit
      EOF

    else
      log :no_results, "No alleles where identified"
    end

    snapshot.glob("*")
  end
end
