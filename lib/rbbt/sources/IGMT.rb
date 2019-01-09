require 'rbbt-util'
require 'rbbt/resource'

module IGMT
  extend Resource
  self.subdir = 'share/databases/IGMT'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  def self.snake_case(allele)
    simple = 'hla_' + allele.downcase.gsub(/[*:]/,'_')
    parts = simple.split("_")
    parts.shift if parts[0] == parts[1]
    parts.last.upcase!
    parts * "_"
  end

  IGMT.claim IGMT["hla_gen.fasta"], :url, "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_gen.fasta"
  IGMT.claim IGMT.identifiers, :proc do
    identifiers = TSV.setup({}, "Allele~Code,Simple,Size#type=:list")
    TSV.traverse IGMT["hla_gen.fasta"].produce, :type => :array do |line|
      next unless line[0] == ">"
      code, name, size = line[1..-1].split(/\s/)
      simple = IGMT.snake_case(name)
      identifiers[name] = [code, simple, size]
    end
    identifiers.to_s
  end
end

iif IGMT["hla_gen.fasta"].produce.find if __FILE__ == $0
iif IGMT.identifiers.produce(true).find if __FILE__ == $0

