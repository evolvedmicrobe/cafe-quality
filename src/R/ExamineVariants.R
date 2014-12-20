dv = read.csv("/Users/nigel/git/cafe-quality/master_full/local_all_variants.csv")
head(dv)

ldv = dv[dv$Ref=="lambda_NEB3011",]
head(ldv)

rep = aggregate(Ref~type+Pos,ldv,length)

head(rep)

ggplot(rep,aes(x=Pos,,y=Ref,colour=type))+geom_point()


aggregate(Pos~type,ldv,function(x) length(x)/nrow(ldv))

gd = d[d$Reference=="lambda_NEB3011" & d$NumPasses>9,]

nd = ldv[(ldv$zmw %in% gd$ZMW) & ldv$type=="Indel",]

nd$homopolymerChar[nd$homopolymerChar=="T"]="A"
nd$homopolymerChar[nd$homopolymerChar=="G"]="C"


res = aggregate(Pos~homopolymerLength+indeltype+homopolymerChar,nd,function(x) length(x)/nrow(nd))
head(res)

res2 = aggregate(Pos~homopolymerLength+homopolymerChar,nd,function(x) length(x)/nrow(nd))
head(res)

pdf("LambdaIndelErrors.pdf",width=6,height=4)
v = ggplot(res,aes(x=homopolymerLength,y=Pos, colour=homopolymerChar))+facet_grid( .~indeltype)+geom_point()+theme_bw(base_size=8)
v = v +labs(x="Homopolymer Length",y="Percentage of Total Errors", title="Indel errors divided by genomic context for reads with >9X Passes")
v = v + scale_colour_discrete(name="Homopolymer Base",labels = c("A/T","G/C"))
v
dev.off()

