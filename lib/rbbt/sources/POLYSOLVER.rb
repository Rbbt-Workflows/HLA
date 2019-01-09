require 'rbbt/util/docker'
module POLYSOLVER

  def self.run(directory, reference = 'b37')
    params = {}

    target_dir = "/job/" << Misc.obj2digest(directory)
    params[:mounts] = {target_dir => directory}

    cmd = "bash /home/polysolver/scripts/shell_call_hla_type #{target_dir}/file.bam Unknown 1 #{reference} STDFQ 0 #{target_dir}"
    Docker.run('sachet/polysolver:v4', cmd, params)
  end
end

