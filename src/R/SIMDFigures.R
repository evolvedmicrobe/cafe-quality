# Ripped off of make recursion figures
library(grid)
library(RColorBrewer)
setwd("/Users/nigel/git/cafe-quality/doc/")
grid.newpage()
pdf_width=7
pdf_height=4
tpl =     "CTGGGTGCGGCGCCATCGCTAC"
read =    "CTGGGTGCGGCGCCATCGCTAC"
rDelTag = "NNNNGNNNNNNGNNNCNNYNCC"
j = nchar(tpl)+1
i = nchar(rDelTag)+1

g_alpha =1
# Make boxes
mkBase <- function(showInit=FALSE, showHEAD=TRUE)  {
  grid.newpage()
  if(showHEAD) {
    
    lay = grid.layout(nrow=4,ncol=1, heights = unit( c(2,1,1,1), c("lines","lines","lines","null"))
      pushViewport(viewport(layout=lay))
      
      
  }
  
  layout=grid.layout(nrow=3,ncol=4, widths=unit(c(5,1,.2,3),c("lines","null","null","mm")),heights=unit(c(3,1,3),c("lines","null","mm")))
  grid.rect()
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
      if(showInit) {
        if (ii==1 | jj==1) {
          if(ii==1 & jj==1)
          {grid.text(0, just=c("centre")) }
          else {
            grid.text(expression(-infinity), just=c("centre"))
          }
        }
      }
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

mkBase(FALSE)


mkLegend <- function(names, colors) {
  seekViewport("base")
  pushViewport(viewport(layout.pos.col=3, layout.pos.row=2, name = "outerlegend"))  
  mat = grid.layout(3, 1)
  pushViewport(viewport(layout=mat))
  pushViewport(viewport(layout.pos.row=2, name = "legendcell"))  
  n = length(names)
  matn = grid.layout(n,4, widths=unit(c(1,.2,.8,1),c("mm","null","null","mm")))
  pushViewport(viewport(layout=matn, name = "legend"))
  for(i in 1:n) {
    seekViewport("legend")
    pushViewport(viewport(layout.pos.row=i, layout.pos.col=2))
    size =unit(.75,"npc")
    grid.rect(width = size, height = size)
    grid.rect(width = size, height = size, gp=gpar(fill=colors[i], alpha = g_alpha))
    upViewport()
    pushViewport(viewport(layout.pos.row=i, layout.pos.col=3))
    grid.text(names[[i]], x=0.05, just=0)    
    upViewport()
  }
  seekViewport("base") 
}


fillRect <-function(i,j,color)
{
  seekViewport("matrix")  
  vp = paste(i,j,sep=",")
  downViewport(vp)
  grid.rect(width=unit(1, "npc"), height=unit(1, "npc"),gp=gpar(fill=color, alpha = 1))
  seekViewport("base")
}

drawTriangle<-function(i,j,color) {
  seekViewport("matrix")  
  vp = paste(i,j,sep=",")
  downViewport(vp)
  x <- c(0, 1, 1)
  y <- c(0, 1, 0)
  grid.path(x,y, gp=gpar(fill=color))
  seekViewport("base")
}

fillDiag <- function(i, j, color, size = 16) {
  els = seq(0, size - 1)
  rows = i + els
  cols = j - els
  for(q in (els+1)) {
    fillRect(rows[q],cols[q], color)
  } 
}

drawArrow<- function(i1,j1,i2,j2) {
  # Grid line to is the most confusing function in R, it 
  # magically appears to remember state, will need to investigate more
  vp1 = paste(i1,j1,sep=",")
  vp2 = paste(i2,j2,sep=",")
  seekViewport(vp1)
  grid.move.to(.5,.5)
  seekViewport(vp2)
  grid.line.to(x=.75,y=.5,arrow=arrow(angle = 30, length = unit(.75, "npc"),ends = "last", type = "open"),gp=gpar(lwd=5))#, vp=current.viewport())
  seekViewport("base")  
}

sr = 5
sc = 19
#pdf("sse_elements.pdf" , width = pdf_width, height = pdf_height)
mkBase(FALSE)
colors = brewer.pal(5,"Set1")# list("green", "red", "blue", "black", "black")
fillDiag(sr,sc, colors[1])
fillDiag(sr-1,sc-1,colors[3])
fillDiag(sr,sc-1,colors[4]) #Deletion
fillDiag(sr-1,sc-2,colors[5])
fillDiag(sr-1,sc, colors[2])
#now to manually add the split
xs = sr + seq(0,14)
ys = sc -1 - seq(0,14)
for (q in seq(1,15)) { drawTriangle(xs[q],ys[q], colors[4])}
names = list("Current", "Insertion", "Match", "Deletion","Merge")
drawArrow(i-1,2,i-1,sc)
drawArrow(sr-3,2,sr-3,sc)
mkLegend( names, colors  )
#dev.off()



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







