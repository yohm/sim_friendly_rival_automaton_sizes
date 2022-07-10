require "fileutils"

if ARGV.size != 2
  puts "Usage: split_lines.rb <file> <num_lines>"
  exit 1
end

num_lines = ARGV[1].to_i

# Split lines in a file into multiple files
dir_idx = 0
file = open(ARGV[0], 'r')

while !file.eof? do
  dir_name = "split_" + sprintf("%03d", dir_idx)
  FileUtils.mkdir_p(dir_name)
  FileUtils.cd(dir_name)
  FileUtils.ln_s("../fugaku_job.sh", ".", :force => true)
  FileUtils.cd("..")

  fout = open("#{dir_name}/#{File.basename(ARGV[0])}", 'w')
  num_lines.times do |i|
    line = file.readline rescue break
    fout.write(line)
  end
  dir_idx += 1
end