require 'rbbt/util/docker'
module POLYSOLVER

  def self.run(directory, reference = 'b37')
    params = {}

    reference = 'hg19' if reference == 'b37'

    target_dir = "/job/" << Misc.obj2digest(directory)
    params[:mounts] = {target_dir => directory}

    cmd = "env SAMTOOLS_DIR=/home/polysolver/binaries/ bash /home/polysolver/scripts/shell_call_hla_type #{target_dir}/file.bam Unknown 1 #{reference} STDFQ 0 #{target_dir}"
    Docker.run('mikisvaz/polysolver:v5', cmd, params)
  end
end

