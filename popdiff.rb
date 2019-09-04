#!/usr/bin/env ruby
require 'optparse'

genotypefile = nil
popfile = nil
missing_rate = 1

opt = OptionParser.new
opt.on('-i Genotype file', 'required') {|v| genotypefile = v}
opt.on('-p Population file', 'required') {|v| popfile = v}
opt.on('-m Missing data allowed in each population', 'default: 1.0') {|v| missing_rate = v.to_f}
opt.on('-h', '--help', 'Show this message') { puts opt ; exit }
opt.parse!(ARGV)

ind2pop = Hash.new
File.open(popfile).each do |line|
    ind = line.slice(/^(.+?)\s.+\n/, 1)
    pop = line.slice(/^.+?\s(.+)\n/, 1)
    next if ind == nil
    ind2pop[ind] = pop
end

num_pop = ind2pop.values.uniq.size
ind = nil
File.open(genotypefile).each do |line|
    a = line.chomp.split("\t")
    a.shift(4)
    if line =~ /^Scaffold_ID/
        puts line.chomp + "\tDifferentiated\tUndifferentiated"
        ind = a
    else
        grouped_gtypes = Hash.new { |h,k| h[k] = [] }
        a.size.times do |i|
            gtype = a[i]
            pop = ind2pop[ind[i]]
            grouped_gtypes[pop] << gtype
        end
        allele_by_pop = [] ; grouped_alleles = []
        Hash[grouped_gtypes.sort].values.each do |x|
            if x.count('-/-')/x.size.to_f > missing_rate
                # puts line.chomp + "\tToo many missing genotype data"
                break
            else
                alleles = x.join('').gsub(/\/|-/, '')
                alleles = alleles.split('').uniq.sort
                if alleles.size > 0
                    ca = alleles.count("A")
                    ct = alleles.count("T")
                    cg = alleles.count("G")
                    cc = alleles.count("C")
                    allele_by_pop << [ca, ct, cg, cc]
                    grouped_alleles << alleles
                end
            end
        end
        next unless allele_by_pop.size == num_pop
        t_ap = allele_by_pop.transpose
        count_pop = t_ap.map {|x| x.inject(:+)}
        next if count_pop.max == num_pop
        # p grouped_alleles
        popcomb = []
        ind2pop.values.uniq.sort.combination(2) do |x,y|
            popcomb << "#{x}-#{y}"
        end
        popdiff = [] ; popundiff = []
        n = 0
        allele_by_pop.combination(2) do |x,y|
            count_pop = [x, y].transpose.map{|z| z.inject(:+)}
            if count_pop.max == 1
                # p grouped_alleles
                popdiff << popcomb[n]
            elsif count_pop.max == 2
                popundiff << popcomb[n]
            end
            n += 1
        end
        if popdiff.size == n
            puts line.chomp + "\t" + "All_populations" + "\t" + "-"
        elsif popdiff.size > 0 && popdiff.size < n
            puts line.chomp + "\t" + popdiff.join(',') + "\t" + popundiff.join(',')
        end
    end
end