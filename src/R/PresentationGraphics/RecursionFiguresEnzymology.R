setwd("/Users/nigel/git/cafe-quality/doc/")
grid.newpage()
pdf_width=6
pdf_height=3


# Make boxes
mkBase <- function(tpl = "CTGGGGA",
                   read =    "CTGGGA",
                   rDelTag = "NNNNGN")  {
  grid.rect(width=unit(1, "npc"), height=unit(1, "npc"),gp=gpar(fill="white"))
  
  j = nchar(tpl)+1
  i = nchar(rDelTag)+1
  layout=grid.layout(nrow=3,ncol=3, widths=unit(c(5,1,3),c("lines","null","mm")),heights=unit(c(3,1,3),c("lines","null","mm")))
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
  #now for top row
  upViewport()
  
  #now for top row
  upViewport()
  
  # Add column bases
  seekViewport("base")
  pushViewport(viewport(layout.pos.col=2, layout.pos.row=1,name="top"))
    tempRow = grid.layout(2, 1)
    pushViewport(viewport(layout=tempRow,name="template"))
    
      pushViewport(viewport(layout.pos.col=1, layout.pos.row=1,name="read"))
      grid.text("Template Sequence",gp=gpar(fontface="bold"))
      upViewport()
      pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
      pushViewport(viewport(layout=grid.layout(1,j),name="TemplateSequence"))
      for(jj in 1:(j-1)) {
        name = paste("tpl",jj,sep="_")
        pushViewport(viewport(layout.pos.col=(jj+1), layout.pos.row=1,name=name))
        grid.text(substr(tpl,jj,jj ))
        upViewport()
      }
      seekViewport("base")
  
  # Add read bases
  seekViewport("base")
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=2,name="readColumns"))  
  tempRow = grid.layout(1, 2)
  pushViewport(viewport(layout=tempRow,name="read_data"))
  
    pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
    pushViewport(viewport(layout=grid.layout(1,2), name="read"))
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
    grid.text("Read Sequence",rot=90,gp=gpar(fontface="bold"))
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
    seekViewport("read_data")
  
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
    pushViewport(viewport(layout=grid.layout(1,2), name="read_del"))
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
    grid.text("Deletion Tags",rot=90,gp=gpar(fontface="bold"))
    upViewport()
    pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
    pushViewport(viewport(layout=grid.layout(i,1),name="DeletionSequence"))
    for(ii in 1:(i-1)) {
      name = paste("del",ii,sep="_")
      pushViewport(viewport(layout.pos.col=1, layout.pos.row=(ii+1),name=name))
      grid.text(substr(rDelTag,ii,ii ))
      upViewport()   
    }
    seekViewport("base")
}



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
    grid.line.to(x=.45,y=.3,arrow=arrow(angle = 30, length = unit(.5, "npc"),ends = "last", type = "open"))#, vp=current.viewport())
  }
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



## Now make the plot
dev.off()
dev.new()
pdf("RecursionMatch.pdf",height=4,width=6)
mkBase()
for(oi in 1:3) {
drawArrow(oi,oi,oi+1,oi+1)
fillRect(oi+1,oi+1,"red")
}
drawArrow(4,4,5,6)
fillRect(5,6,"blue")
for(oi in 5:6) {
  drawArrow(oi,oi+1,oi+1,oi+2)
  fillRect(oi+1,oi+2,"red")
}
dev.off()
