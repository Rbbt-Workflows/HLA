require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

require 'rbbt/sources/IGMT'
require 'rbbt/sources/SOAPHLA'
require 'rbbt/sources/opti_type'
require 'rbbt/sources/LOHHLA'
require 'rbbt/sources/POLYSOLVER'
require 'rbbt/sources/xHLA'

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

#require 'rbbt/knowledge_base/HLA'
#require 'rbbt/entity/HLA'

require 'HLA/tasks/sample' if defined? Sample
