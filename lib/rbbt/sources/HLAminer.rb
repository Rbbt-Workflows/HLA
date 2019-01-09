require 'rbbt-util'
require 'rbbt/resource'

module HLAminer
  extend Resource
  self.subdir = 'share/databases/HLAminer'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib


  Rbbt.claim Rbbt.software.opt.HLAminer, :install, Rbbt.share.install.software.HLAminer.find
end

iif Rbbt.software.opt.HLAminer.find if __FILE__ == $0
