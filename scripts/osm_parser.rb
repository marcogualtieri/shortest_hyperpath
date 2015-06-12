require 'OSM/objects'
require 'OSM/StreamParser'
require_relative 'wgs84.rb'

$node_counter = 0
$way_counter = 0
$nodes = Array.new
$ways = Array.new
$nodes_hash = Hash.new

class MyCallbacks<OSM::Callbacks
  def node(node)
    stop = 0
    name = node.id.to_s
    $nodes[$node_counter] = $node_counter.to_s + ' ' + stop.to_s + ' ' + node.lat.to_s + ' ' + node.lon.to_s + ' ' + name
    h = {:id => $node_counter, :lat => node.lat, :lon => node.lon}
    $nodes_hash[node.id] = h
    $node_counter = $node_counter + 1
  end

  def way(way)
    #name  
    name = 'unnamed'
    if(!way.name.nil?) then name = way.name end
    #nodes
    for index in 0..(way.nodes.length-2)
      i = $nodes_hash[way.nodes[index].to_i]
      j = $nodes_hash[way.nodes[index+1].to_i]
      if(i.nil? or j.nil?)
        puts "Error in way #{way.id} - #{way.name}: node(s) not found ( #{way.nodes[index]} , #{way.nodes[index+1]} )"
      else
        meters = vicenti_distance(i[:lat].to_f,i[:lon].to_f,j[:lat].to_f,j[:lon].to_f)
        weight = (meters / 1.5).round(3) #m/s: velocità a passo d'uomo (5-6 km/h) 
        $ways[$way_counter] = i[:id].to_s + ' ' + j[:id].to_s + ' ' + weight.to_s + ' ' + name
        $way_counter = $way_counter + 1
        $ways[$way_counter] = j[:id].to_s + ' ' + i[:id].to_s + ' ' + weight.to_s + ' ' + name
        $way_counter = $way_counter + 1
      end
    end      
  end
end

#main
if ARGV.length != 1 
  puts "usage: $ruby osm_parser.rb filename.osm";
  exit
end
filename = ARGV[0]
new_filename = filename.split('.')[0] + '_parsed'
puts 'Parsing ' + filename + ' to ' + new_filename + '...'

cb = MyCallbacks.new
parser = OSM::StreamParser.new(filename: filename, callbacks: cb)
parser.parse

File.open(new_filename,'w') do |f|
  f.puts $node_counter
  $nodes.each do |line|
    f.puts line
  end
  f.puts
  f.puts $way_counter
  $ways.each do |line|
    f.puts line
  end
end
