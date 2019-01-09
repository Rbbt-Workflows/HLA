require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/IGMT'
require 'rbbt/sources/SOAPHLA'
require 'rbbt/sources/opti_type'
require 'rbbt/sources/LOHHLA'
require 'rbbt/sources/POLYSOLVER'

Workflow.require_workflow "HTS"
require 'tools/IGV'

module HLA
  extend Workflow

  HLA_REFERENCE = IGMT["hla_gen.fasta"].produce.find

  Rbbt::Config.set(:cpus, 3, "razers::20")

end

require 'HLA/tasks/reads.rb'
require 'HLA/tasks/typing.rb'
require 'HLA/tasks/alignment.rb'

Workflow.require_workflow "Sample"
Workflow.require_workflow "PVacSeq"

module Sample

  dep :BAM
  dep :reference
  dep HLA, :SOAPHLA, :BAM => :BAM, :reference => :reference
  task :hla_genotype => :array do
    Open.read(step(:SOAPHLA).path).split("\n").collect{|line|
      next if line == '---'
      line.split("\t").values_at 0, 1
    }.flatten.compact
  end

  dep :genomic_mutations
  dep :hla_genotype
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
#require 'rbbt/knowledge_base/HLA'
#require 'rbbt/entity/HLA'

