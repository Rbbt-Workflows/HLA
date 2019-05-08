Workflow.require_workflow "PVacSeq"

module Sample
  dep :BAM
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg19 hg38 GRCh38 hs37d5), :nofile => true
  dep HLA, :SOAPHLA, :BAM => :BAM, :reference => :reference
  task :hla_genotype => :array do
    Open.read(step(:SOAPHLA).path).split("\n").collect{|line|
      next if line == '---'
      line.split("\t").values_at 0, 1
    }.flatten.compact
  end

  dep :genomic_mutations, :compute => :produce
  dep :hla_genotype, :compute => :produce
  dep PVacSeq, :analysis, :positions => :genomic_mutations, :alleles => :hla_genotype
  task :neo_epitopes => :tsv do

    parser = TSV::Parser.new step(:analysis).join, :type => :list
    dumper = TSV::Dumper.new parser.options.merge(:key_field => "Genomic Mutation", :namespace => CFDTL.organism, :fields => parser.fields[4..-1])
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
