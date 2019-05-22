require 'rbbt/util/docker'
module XHLA

  def self.run(directory, reference = 'b37')
    params = {}

    reference = 'hg19' if reference == 'b37'

    target_dir = "/job/" << Misc.obj2digest(directory)
    params[:mounts] = {target_dir => directory}

    cmd = "    --sample_id sample --input_bam_path '#{target_dir}/file.bam' \
        --output_path '#{target_dir}'"
    Docker.run('humanlongevity/hla', cmd, params)
  end
end

