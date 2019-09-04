dir = ARGV[0]

tag_cat_file = Dir.glob("#{dir}/batch_*.catalog.tags.tsv")[0]
snp_cat_file = Dir.glob("#{dir}/batch_*.catalog.snps.tsv")[0]
tagfiles = Dir.glob("#{dir}/[^batch]*.tags.tsv").sort
snpfiles = Dir.glob("#{dir}/[^batch]*.snps.tsv").sort
mtcfiles = Dir.glob("#{dir}/[^batch]*.matches.tsv").sort

dataframe = Hash.new

sample_size = snpfiles.size

fname2sid = {}
snpfiles.each do |file|
    File.open(file) do |f|
        second_line = ""
        2.times{second_line = f.gets}
        next unless second_line
        sid = second_line.split("\t")[1]
        fname2sid[file] = sid
    end
end

taginfo_cat = {}
tag_sam2cat = {}
File.open(tag_cat_file).each do |line|
    next if line =~ /^#/
    a = line.split("\t")
    tagcat = a[2]
    scaf = "-"
    start = a[4].to_i
    strand = a[5]
    taginfo_cat[tagcat] = [scaf, start, strand]
    taglink = a[8]
    taglink.split(',').each {|l| tag_sam2cat[l] = tagcat}
end

def taginfo(tagfile)
    hash = {}
    File.open(tagfile).each do |line|
        next unless line =~ /consensus/
        a = line.split("\t")
        sampleID = a[1]
        tag = a[2]
        scaf = "-"
        start = a[4].to_i
        strand = a[5]
        id = sampleID + "_" + tag
        hash[id] = [scaf, start, strand]
    end
    return hash
end

def comp_base(ary)
    a = ary.join.upcase.gsub(/A/, 't')
    b = a.gsub(/T/, 'A')
    c = b.gsub(/C/, 'g')
    d = c.gsub(/G/, 'C').upcase
    return d.split('')
end

File.open(snp_cat_file).each do |line|
    next if line =~ /^#/
    a = line.split("\t")
    tag = a[2]
    column = a[3].to_i
    info = taginfo_cat[tag]
    scaf = info[0]
    start = info[1].to_i
    strand = info[2]
    pos = strand == '+' ? start+1+column : start+1-column
    scaf_pos_tag = scaf + "_" + pos.to_s + "_" + tag
    dataframe[scaf_pos_tag] = [scaf, pos, tag] + ["-/-"] * sample_size
end

header = "Scaffold_ID\tPosition\tTag_ID\tN_Alleles"

sample_size.times do |i|

    sample_name = tagfiles[i].slice(/.+?\/(.+?)\.tags\.tsv/, 1)
    header << "\t" + sample_name
    taginfo_sample = taginfo(tagfiles[i])

    File.open(snpfiles[i]).each do |line|
        next if line =~ /^#/
        a = line.split("\t")
        tagsam = a[2]
        sampleID = fname2sid[snpfiles[i]]
        samtagid = sampleID + "_" + tagsam
        tagcat = tag_sam2cat[samtagid]
        next unless tagcat
        column = a[3].to_i
        info = taginfo_sample[samtagid]
        scaf = info[0]
        start = info[1]
        strand = info[2]
        pos = strand == '+' ? start+1+column : start+1-column
        scaf_pos_tag = scaf + "_" + pos.to_s + "_" + tagcat
        data = dataframe[scaf_pos_tag]
        if data
            type = a[4]
            allele1 = a[6]
            allele2 = a[7]
            genotype =
                case type
                when "O" ; [allele1, allele1]
                when "E" ; [allele1, allele2]
                when "U" ; ["-", "-"]
                else ; ["-", "-"]
                end
            genotype = comp_base(genotype) if strand == '-'
            data[i+3] = genotype.sort.insert(1, '/').join
            dataframe[scaf_pos_tag] = data
        end
    end

end

puts header

dataframe.sort.each do |k, v|
    vd = v.dup
    vd.shift(3)
    vd.delete("-/-")
    alleles = []
    vd.each {|g| alleles << g.split('/')}
    num_alleles = alleles.flatten.uniq.size
    if num_alleles > 1
        v.insert(3, num_alleles)
        puts v.join("\t")
    end
end