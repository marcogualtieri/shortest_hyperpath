#!/usr/bin/ruby
require 'rubygems'
require 'shapelib'

if ARGV.length != 1 
  puts "usage: $ruby hypergraph_to_shape.rb filename";
  exit
end

puts "Opening " + ARGV[0]

nodes = {}

sfe = Shapelib::ShapeFile::new(ARGV[0] + "_edge.shp", :Arc, 
                               [
                                ['i', :Integer, 12],
                                ['j', :Integer, 12],
                                ['cost', :String, 64],
				['prob', :String, 64],
				['name', :String, 64]
                               ])

sfv = Shapelib::ShapeFile::new(ARGV[0] + "_vertex.shp", :Point, 
                               [
				['id', :Integer, 12],
				['stop_vertex', :Integer, 12],
				['exp_cost', :String, 64],
				['comb_freq', :String, 64],
				['name', :String, 64]
			       ])

file = File.open(ARGV[0])
lines = file.readlines
#vertex
n = lines[0].split[0].to_i
for index in 1..n
  elem = lines[index].split
  nodes[elem[0]] = [elem[2],elem[3]]
  name = ''
  for k in 6..elem.size-1
    name = name + elem[k] + ' ' 
  end
  v = Shapelib::Point::new(elem[3].to_f,
                             elem[2].to_f)
  v["id"] = elem[0].to_i
  v["stop_vertex"] = elem[1].to_i
  v["exp_cost"] = elem[4]
  v["comb_freq"] = elem[5]
  v["name"] = name
  sfv.write(v)
end
#edge
m = lines[n+2].split[0].to_i
for index in n+3..n+2+m
  elem = lines[index].split
  i = nodes[elem[0]]
  j = nodes[elem[1]]
  if(i.nil?)
    puts "Node " + elem[0] + " not found!"
  elsif(j.nil?)
    puts "Node " + elem[1] + " not found!"
  else
    name = ''
    for k in 4..elem.size-1
      name = name + elem[k] + ' ' 
    end
    e = Shapelib::Arc::new(
           Shapelib::Point::new(i[1].to_f, 
                                i[0].to_f),
           Shapelib::Point::new(j[1].to_f, 
                                j[0].to_f))
    e["i"] = elem[0].to_i
    e["j"] = elem[1].to_i
    e["cost"] = elem[2]
    e["prob"] = elem[3]
    e["name"] = name
    sfe.write(e)
  end
end
