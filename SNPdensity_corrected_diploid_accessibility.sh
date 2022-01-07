samtools view -@ ${threads} ${junction_rt_string} -h reads.sorted.bam $chr | awk -v chr=$chr 'BEGIN{
        OFS="\t"}FILENAME==ARGV[1]{
                if($2==chr"-r"||$2==chr"-a"){
                        if(keep[$1]&&keep[$1]!=$2){
                                delete keep[$1]
                        }else{
                                if (keep[$1]){
                                        x=rand();
                                        if (x>=0.5){
                                                keep[$1]=$2;
                                                keepRT[$1]=$3;
                                        }
                                } else {
                                        keep[$1]=$2;
                                        if ($3%2==0){
                                                keepRT[$1]=$3+1;
                                        } else{
                                                keepRT[$1]=$3-1;
                                        }
                                }
                        }
                };
                next
        }
        $0~/^@/{next}
        ($1 in keep){
                $3=keep[$1]; 
                
                for (i=12; i<=NF; i++) {
                        if ($i~/^ip/) {
                                split($i, ip, ":");
                        }
                        else if ($i ~ /^rt:/) {
                                split($i, rt, ":");
                        }
                }
                if (rt[3]==keepRT[$1]) {
                        locus[$3" "ip[3]]++;
                }
        }END{
                for (i in locus) {
                        split(i, a, " "); 
                        print a[1], a[2]-1, a[2], locus[i]
                }
        }' reads_to_homologs.txt - | sort -k1,1 -k2,2n --parallel=${threads} -S 6G > tmp.bedgraph
