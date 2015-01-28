setwd("/Users/nigel/git/cafe-quality/doc/PBEP_3/")
library(grid)

grid.newpage()
pdf_width=6
pdf_height=3
tpl =    "CTGGGGA"
read =    "CTGGGA"
rDelTag = "NNNNGN"
j = nchar(tpl)+1
i = nchar(rDelTag)+1

# Make boxes
mkBase <- function()  {
  layout=grid.layout(nrow=3,ncol=3, widths=unit(c(5,1,3),c("lines","null","mm")),heights=unit(c(3,1,3),c("lines","null","mm")))
  # Layout is
  # 1-1
  # 1-2 Template sequence
  # 2-1 Read sequence
  # 2-2 matrix
  pushViewport(viewport(layout=layout,name="base"))
  pushViewport(viewport(layout.pos.col=2, layout.pos.row=2, name = "matrix"))  
  mat = grid.layout(i, j)
  #grid.show.layout(mat)
  pushViewport(viewport(layout=mat,name="mat"))
  grid.rect()
  #grid.grill()
  for (ii in 1:i) {
    for(jj in 1:j) {
      name=paste(ii,jj,sep=",")    
      pushViewport(viewport(layout.pos.col=jj, layout.pos.row=ii, name = name)) 
      grid.rect()
      upViewport()
    }
  }

  
  # Add column bases
  seekViewport("base")
  pushViewport(viewport(layout.pos.col=2, layout.pos.row=1,name="top"))
  tempRow = grid.layout(2, 1)
  pushViewport(viewport(layout=tempRow,name="read"))
  
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=1,name="read"))
  grid.text("Read Sequence",gp=gpar(fontface="bold"))
  upViewport()
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
  pushViewport(viewport(layout=grid.layout(1,j),name="ReadSequence"))
  for(jj in 1:(j-1)) {
    name = paste("tpl",jj,sep="_")
    pushViewport(viewport(layout.pos.col=(jj+1), layout.pos.row=1,name=name))
    grid.text(substr(tpl,jj,jj ))
    upViewport()
  }
  
  # Add row bases
  seekViewport("base")
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=2,name="readColumns"))  
  tempRow = grid.layout(1, 1)
  pushViewport(viewport(layout=tempRow,name="read_data"))
  
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
  pushViewport(viewport(layout=grid.layout(1,2), name="read"))
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
  grid.text("Template Sequence",rot=90,gp=gpar(fontface="bold"))
  upViewport()
  pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
  pushViewport(viewport(layout=grid.layout(i,1),name="ReadSequence"))
  for(ii in 1:(i-1)) {
    ypos = (ii-1)*1/i
    name = paste("read",ii,sep="_")
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=(ii+1),name = name))
    grid.text(substr(read,ii,ii ))
    upViewport()   
  }

    seekViewport("base")
}

mkBase(TRUE)

fillRect <-function(i,j,color)
{
  seekViewport("matrix")  
  vp = paste(i,j,sep=",")
  downViewport(vp)
  grid.rect(width=unit(1, "npc"), height=unit(1, "npc"),gp=gpar(fill=color,alpha=.3))
  seekViewport("base")
}
drawArrow<- function(i1,j1,i2,j2,noise=FALSE) {
  # Grid line to is the most confusing function in R, it 
  # magically appears to remember state, will need to investigate more
  vp1 = paste(i1,j1,sep=",")
  vp2 = paste(i2,j2,sep=",")
  seekViewport(vp1)
  grid.move.to(.5,.3)
  seekViewport(vp2)
  if(noise) {
    grid.line.to(x=rnorm(.5,.1),y=.3,arrow=arrow(angle = 30, length = unit(, "npc"),ends = "last", type = "open"))#, vp=current.viewport())
  }else{
    grid.line.to(x=.5,y=.3,arrow=arrow(angle = 30, length = unit(.5, "npc"),ends = "last", type = "open"))#, vp=current.viewport())
  }
  seekViewport("base")  
}
drawTextArrow<- function(box1,box2) {
  # Grid line to is the most confusing function in R, it 
  # magically appears to remember state, will need to investigate more
  seekViewport(box1)
  grid.move.to(.5,.7)
  seekViewport(box2)
  grid.line.to(x=.4,y=.1,arrow=arrow(angle = 30, length = unit(.5, "npc"),ends = "both", type = "open"),gp=gpar(lwd=2))#, vp=current.viewport())
  seekViewport("base")  
}
addText <- function(ii,jj,text) {
  vp1 = paste(ii,jj,sep=",")
  seekViewport(vp1)
  grid.text(text)
  seekViewport("base") 
}
colorLetter <- function(ind, prefix, color="blue") {
  # Acceptable prefixes include "read_", del_, tpl_
  seq="N"
  if (prefix=="read") {seq=read} else if (prefix=="del") {seq = rDelTag} else if (prefix=="tpl") {seq=tpl}
  name = paste(prefix,ind,sep="_")
  seekViewport(name)
  grid.text(substr(seq,ind,ind ), gp=gpar(col=color))
  seekViewport("base")
}

pdf("Initialization.pdf", width=pdf_width, height= pdf_height)
mkBase(TRUE)
dev.off()

## Now make the insertion move plot
pdf("Insertion.pdf", width=pdf_width, height= pdf_height)
mkBase()
oi = 5
oj = 5
fillRect(oi,oj,"blue")
fillRect(oi-1,oj,"red")
addText(oi-1,oj,expression(alpha["i-1,j"]))
drawArrow(oi-1,oj,oi,oj)
colorLetter(oi-1,"read")
colorLetter(oj,"tpl","red")
drawTextArrow("read_4","tpl_5")
dev.off()

## Now make the deletion move plot
pdf("Deletion.pdf", width=pdf_width, height= pdf_height)
grid.newpage()
mkBase()
fillRect(oi,oj,"blue")
fillRect(oi,oj-1,"red")
addText(oi,oj-1,expression(alpha["i,j-1"]))
drawArrow(oi,oj-1,oi,oj)
colorLetter(oj-1,"tpl","red")
colorLetter(oj,"tpl","red")
drawTextArrow("tpl_4","tpl_5")
dev.off()

#Now the Match plot
pdf("Match.pdf", width=pdf_width, height= pdf_height)
grid.newpage()
mkBase()
fillRect(oi,oj,"blue")
fillRect(oi-1,oj-1,"red")
addText(oi-1,oj-1,expression(alpha["i-1,j-1"]))
drawArrow(oi-1,oj-1,oi,oj)
colorLetter(oi-1,"read")
colorLetter(oj-1,"tpl","red")
drawTextArrow("read_4","tpl_4")
dev.off()






