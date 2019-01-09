require 'rbbt-util'
require 'rbbt/resource'
require_relative 'IGMT'

module LOHHLA
  
  extend Resource
  self.subdir = 'share/databases/LOHHLA'

  Rbbt.claim Rbbt.software.opt.Jellyfish, :install, Rbbt.share.install.software.Jellyfish.find
  Rbbt.claim Rbbt.software.opt.LOHHLA, :install, Rbbt.share.install.software.LOHHLA.find
  Rbbt.claim Rbbt.software.opt.PicardTools, :install, Rbbt.share.install.software.PicardTools.find

  LOHHLA.claim LOHHLA.fasta["fixed_igmt.fasta"], :proc do
    code2simple = IGMT.identifiers.index :target => "Simple"
    TSV.traverse IGMT["hla_gen.fasta"].produce.find, :into => :stream, :type => :array do |line|
      if line[0] == '>'
        code = line.split(/\s/).first[1..-1]
        simple = code2simple[code] || code
        ">" << simple
      else
        line
      end
    end
  end

  def self.run(output, alleles, normalBAM, otherBAMs, reference = 'b37', purity = 1.0, ploidy = 2.0)

    output = output.find if Path === output
    bam_dir = File.join(output, 'data', 'bam')
    hla_path = File.join(output, 'data', 'hlas')
    copy_number = File.join(output, 'data', 'cnv')

    Misc.in_dir bam_dir do
      CMD.cmd("ln -s '#{normalBAM}'* .", :no_fail => true)
      otherBAMs.each do |other|
        CMD.cmd("ln -s '#{other}'* .", :no_fail => true)
      end
    end

    alleles = alleles.collect{|a| IGMT.snake_case(a) }

    Open.write(hla_path, alleles * "\n" + "\n")

    sample = otherBAMs.collect{|other| File.basename(other) }.first.sub('.bam','')

    Open.write(copy_number, <<-EOF)
Ploidy	tumorPurity	tumorPloidy
#{sample}	2	#{purity}	#{ploidy}
    EOF

    normal_bam_local = File.join(bam_dir, File.basename(normalBAM))

    fasta_loc = Samtools.prepare_FASTA(LOHHLA.fasta["fixed_igmt.fasta"])

    lohhla_script = if reference == 'b37'
                      Rbbt.software.opt.LOHHLA.produce["LOHHLAscript.b37.R"].find
                    else
                      Rbbt.software.opt.LOHHLA.produce["LOHHLAscript.R"].find
                    end

    script = <<-EOF

alias jellyfish='#{Rbbt.software.opt.bin.jellyfish}';
#alias samtools='#{Samtools::Samtools_CMD}';
#alias bedtools='#{Rbbt.software.opt.bin.bedtools}';

Rscript #{lohhla_script} \
  --gatkDir '#{Rbbt.software.opt.PicardTools.find}' --novoDir '#{Rbbt.software.opt.NovoAlign.find}' \
  --outputDir '#{output}' \
  --BAMDir '#{bam_dir}' \
  --normalBAMfile '#{normal_bam_local}' \
  --patientId "Sample" \
  --hlaPath '#{hla_path}' \
  --HLAfastaLoc '#{fasta_loc}' \
  --HLAexonLoc '#{Rbbt.software.opt.LOHHLA.produce.data["hla.dat"].find}' \
  --CopyNumLoc #{copy_number} --mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE  
    EOF

    CMD.cmd_log(script)
  end
end

if __FILE__ == $0
  Log.severity = 0
  iif Rbbt.software.opt.NovoAlign.produce
  iif Rbbt.software.opt.Jellyfish.produce
  iif Rbbt.software.opt.LOHHLA.produce
  iif Rbbt.software.opt.PicardTools.produce
  iif LOHHLA.fasta["fixed_igmt.fasta"].produce
end
