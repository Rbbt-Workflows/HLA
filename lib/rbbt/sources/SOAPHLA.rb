require 'rbbt-util'
require 'rbbt/resource'

module SOAPHLA
  extend Resource
  self.subdir = 'share/databases/SOAPHLA'

  #def self.organism(org="Hsa")
  #  Organism.default_code(org)
  #end

  #self.search_paths = {}
  #self.search_paths[:default] = :lib

  def self.run(file, dir, build = 'b37')
    version = 'hg19' if build == 'b37'
    tmp_file = file
    Misc.in_dir Rbbt.software.opt.SOAPHLA.produce.find do 
      name = File.basename(tmp_file, '.bam')
      begin
        if build == 'b37'
          CMD.cmd_log("perl MHC_autopipeline.b37.pl -i '#{tmp_file}' -od '#{dir}' -v #{version}")
        else
          CMD.cmd_log("perl MHC_autopipeline.pl -i '#{tmp_file}' -od '#{dir}' -v #{version}")
        end
      rescue ProcessFailed
        raise RbbtException, "Error running SOAPHLA command"
      end
    end
  end

  Rbbt.claim Rbbt.software.opt.SOAPHLA, :install, Rbbt.share.install.software.SOAPHLA.find

end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.SOAPHLA.produce 
end
