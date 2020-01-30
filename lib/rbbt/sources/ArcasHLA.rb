require 'rbbt/util/docker'

module ArcasHLA
  Rbbt.claim Rbbt.software.opt.Kallisto, :install, "https://github.com/pachterlab/kallisto.git"
  Rbbt.claim Rbbt.software.opt.ArcasHLA, :install do
    commands =<<-EOF
echo Getting reference
(cd $(opt_dir $name); ./arcasHLA reference  --version 3.34.0)
echo Reference prepared
    EOF
    {:git => "https://github.com/RabadanLab/arcasHLA.git", :commands => commands, :name => 'ArcasHLA'}
  end

  def self.run(bam, alleles, output, cpus = nil)
    cpus ||= Rbbt::Config.get("cpus", :ArcasHLA, :arcashla, :arcas, :default => 1)

    Rbbt.software.opt.Kallisto.produce
    Rbbt.software.opt.ArcasHLA.produce

    Open.mkdir(output)
    name = File.basename(bam).sub('.bam','')
    Misc.in_dir(Rbbt.software.opt.ArcasHLA.find) do
      CMD.cmd_log("./arcasHLA extract '#{bam}' -o '#{output}' --paired -t #{cpus} -v")
      CMD.cmd_log("./arcasHLA genotype #{output}/#{name}.extracted.1.fq.gz #{output}/#{name}.extracted.2.fq.gz -g #{alleles*","} -o #{output} -t #{cpus} -v")
    end
    JSON.parse(Open.read(Dir.glob("#{output}/*.genotype.json").first))
  end
  
  def self.run_fasta(fasta1, fasta2, alleles, output, cpus = nil)
    cpus ||= Rbbt::Config.get("cpus", :ArcasHLA, :arcashla, :arcas, :default => 1)

    Rbbt.software.opt.Kallisto.produce
    Rbbt.software.opt.ArcasHLA.produce

    Open.mkdir(output)
    Misc.in_dir(Rbbt.software.opt.ArcasHLA.find) do
      CMD.cmd_log("./arcasHLA genotype #{fasta1} #{fasta2} -g #{alleles*","} -o #{output} -t #{cpus} -v")
    end
    JSON.parse(Open.read(Dir.glob("#{output}/*.genotype.json").first))
  end
end

