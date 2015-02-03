setwd("/Users/nigel/git/cafe-quality/doc/")
library(grid)
library(gridExtra)
library(wesanderson)

pdf("Hmm.pdf", width=6.5, height=5)
grid.newpage()
tpl =    "AGG"
j = nchar(tpl)

colors =  wes_palette(name = "Moonrise2", type="discrete")
# Layout Branches
cvp = viewport(width=0.95, height=0.95, name="base")
pushViewport(cvp)
grid.rect()

getPoses <- function(j) {
  w = 1 / (2*j)
  s = seq(0, 1, 1/j)
  xs = w + s
  return(xs[1:(length(xs)-1)])
}
xs = getPoses(j) # column centers
ys = rev(getPoses(6)) # row centers

# ADD TEXT ON TOP
for(i in 1:j)
{
  #pushViewport(viewport(layout.pos.col=i, layout.pos.row=1)) 
  txt = substr(tpl,i,(i))
  grid.text(txt,x = unit(xs[i],"npc"), y = unit(ys[1],"npc"), just=c("centre"), gp=gpar(fontsize=25, fontstyle="bold"))
}

#Now add circles for states
labels = c("Stick", "Branch", "Match", "Dark", "Merge")
colors2 = c(colors[1], colors[1], colors[2], colors[3], colors[3])
r = diff(ys)[1]/3
for (ii in 2:6) {
  for(jj in 1:j) {
    name=paste(ii,jj,sep="-")   
    print(name)    
    grid.circle(y = unit(ys[ii], "npc") ,x = unit(xs[jj],"npc"), r=unit(r, "npc"), gp=gpar(fill=colors2[ii-1], alpha=.75), name = name)
    grid.text(y = unit(ys[ii], "npc") ,x = unit(xs[jj],"npc"),labels[ii-1], gp = gpar(fontsize=10))        
  }
}



# NOW ADD ALL POSSIBLE TRANSITIONS
drawTransitions <- function(start_col) {
  
  e_col = start_col+1
  s_col = start_col
  arr <- arrow(angle=25, length=unit(4, "mm"))
  lwd = 2
  
  min_weight = .3 # Amount by which to decrease the thickness of rare events
  
  # STICK TO MATCH
  s_c = paste("2",s_col, sep="-")
  e_c = paste("4", e_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 330), grobY(grid.get(s_c), 330),
                grobX(grid.get(e_c), 120), grobY(grid.get(e_c), 120),
                angle=0, shape= 0, square=TRUE, curvature=0, open=TRUE, ncp=1,
                arrow=arr, gp=gpar(col=colors2[3], lwd= lwd))
  
  # Stick to Merge
  s_c = paste("2",s_col, sep="-")
  e_c = paste("6", e_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 275), grobY(grid.get(s_c), 275),
             grobX(grid.get(e_c), 110), grobY(grid.get(e_c), 120),
             angle=0, shape= 0, square=TRUE, curvature=0, open=TRUE, ncp=1,
             arrow=arr, gp=gpar(col=colors2[4], lwd= lwd * min_weight))
  
  
  # BRANCH TO MATCH
  s_c = paste("3",s_col, sep="-")
  e_c = paste("4",e_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 0), grobY(grid.get(s_c), 0),
             grobX(grid.get(e_c), 150), grobY(grid.get(e_c), 150),
             angle=0, shape= 0, square=TRUE, curvature=0, open=TRUE, ncp=1,
             arrow=arr, gp=gpar(col=colors2[3], lwd=lwd))
  
  # MATCH TO MATCH
  s_c = paste("4",s_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 0), grobY(grid.get(s_c), 0),
             grobX(grid.get(e_c), 180), grobY(grid.get(e_c), 180),
             angle=0, shape= 0, square=TRUE, curvature=0, open=TRUE, ncp=1,
             arrow=arr, gp=gpar(col=colors2[3], lwd=lwd))
  
  # dark to match
  s_c = paste("5",s_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 0), grobY(grid.get(s_c), 0),
             grobX(grid.get(e_c), 210), grobY(grid.get(e_c), 210),
             angle=0, shape= 0, square=TRUE, curvature=0, open=TRUE, ncp=1,
             arrow=arr, gp=gpar(col=colors2[3], lwd=lwd))
  
  # merge to match
  s_c= s_c = paste("6",s_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 0), grobY(grid.get(s_c), 0),
             grobX(grid.get(e_c), 240), grobY(grid.get(e_c), 240),
             angle=0, shape= 0, square=TRUE, curvature=0, open=TRUE, ncp=1,
             arrow=arr, gp=gpar(col=colors2[3], lwd=lwd))
  
  # Match to dark
  s_c = paste("4",s_col, sep="-")
  e_c = paste("5", e_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 330), grobY(grid.get(s_c), 330),
             grobX(grid.get(e_c), 180), grobY(grid.get(e_c), 180),
             angle=0, shape= 0, square=TRUE, curvature=0, open=TRUE, ncp=1,
             arrow=arr, gp=gpar(col=colors2[4], lwd=lwd))
  
  # Match to merge
  e_c = paste("6", e_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 300), grobY(grid.get(s_c), 300),
             grobX(grid.get(e_c), 150), grobY(grid.get(e_c), 150),
             angle=0, shape= 0, square=TRUE, curvature=0, open=TRUE, ncp=1,
             arrow=arr, gp=gpar(col=colors2[4], lwd=lwd))
  
  # Stick to stick
  s_c = paste("2",s_col, sep="-")
  e_c = s_c
  #note decrease curvature for a bigger loop, best to keep square=FALSE to avoid a choppy look
  grid.curve(grobX(grid.get(s_c), 180), grobY(grid.get(s_c), 180),
             grobX(grid.get(e_c), 120), grobY(grid.get(e_c), 120),
             angle=90, shape= 1, square=FALSE,  curvature=-2, open=TRUE, ncp=6,
             arrow=arr, gp=gpar(col=colors2[1], lwd=lwd * min_weight))
  
  # Branch to branch
  s_c = paste("3",s_col, sep="-")
  e_c = s_c
  grid.curve(grobX(grid.get(s_c), 180), grobY(grid.get(s_c), 180),
             grobX(grid.get(e_c), 120), grobY(grid.get(e_c), 120),
             angle=90, shape= 1, square=FALSE,  curvature=-2, open=TRUE, ncp=6,
             arrow=arr, gp=gpar(col=colors2[1], lwd=lwd * min_weight))
  
  # match to stick
  s_c = paste("4",s_col, sep="-")
  e_c = paste("2", s_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 20), grobY(grid.get(s_c), 20),
             grobX(grid.get(e_c), 315), grobY(grid.get(e_c), 315),
             angle=60, shape= 1, square=FALSE,  curvature=.5, open=TRUE, ncp=6,
             arrow=arr, gp=gpar(col=colors2[1], lwd=lwd))
  
  # match to branch
  s_c = paste("4",s_col, sep="-")
  e_c = paste("3", s_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 45), grobY(grid.get(s_c), 45),
             grobX(grid.get(e_c), 290), grobY(grid.get(e_c), 290),
             angle=60, shape= 1, square=FALSE,  curvature=.5, open=TRUE, ncp=6,
             arrow=arr, gp=gpar(col=colors2[1], lwd=lwd))
  
  # Dark to stick
  s_c = paste("5",s_col, sep="-")
  e_c = paste("2", s_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 20), grobY(grid.get(s_c), 20),
             grobX(grid.get(e_c), 30), grobY(grid.get(e_c), 30),
             angle=90, shape= 1, square=FALSE,  curvature=.65, open=TRUE, ncp=8,
             arrow=arr, gp=gpar(col=colors2[1], lwd=lwd * min_weight))
  
  
  # Dark to Branch
  s_c = paste("5",s_col, sep="-")
  e_c = paste("3", s_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 45), grobY(grid.get(s_c), 45),
             grobX(grid.get(e_c), 330), grobY(grid.get(e_c), 330),
             angle=90, shape= 1, square=FALSE,  curvature=.65, open=TRUE, ncp=8,
             arrow=arr, gp=gpar(col=colors2[1], lwd=lwd * min_weight))
  
  # Merge to stick
  s_c = paste("6",s_col, sep="-")
  e_c = paste("2", s_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 20), grobY(grid.get(s_c), 20),
             grobX(grid.get(e_c), 0), grobY(grid.get(e_c), 0),
             angle=90, shape= 1, square=FALSE,  curvature=.45, open=TRUE, ncp=8,
             arrow=arr, gp=gpar(col=colors2[1], lwd=lwd * min_weight))
  
  # Stick to dark
  s_c = paste("2",s_col, sep="-")
  e_c = paste("5", e_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 290), grobY(grid.get(s_c), 290),
             grobX(grid.get(e_c), 130), grobY(grid.get(e_c), 130),
             angle=90, shape= 1, square=FALSE,  curvature=0, open=TRUE, ncp=1,
             arrow=arr, gp=gpar(col=colors2[4], lwd=lwd * min_weight))
  
  # dark to merge
  s_c = paste("5",s_col, sep="-")
  e_c = paste("6", e_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 330), grobY(grid.get(s_c), 330),
             grobX(grid.get(e_c), 180), grobY(grid.get(e_c), 180),
             angle=90, shape= 1, square=FALSE,  curvature=0, open=TRUE, ncp=1,
             arrow=arr, gp=gpar(col=colors2[4], lwd=lwd * min_weight))
  
  # merge to merge
  s_c = paste("6",s_col, sep="-")
  e_c = paste("6", e_col, sep="-")
  grid.curve(grobX(grid.get(s_c), 0), grobY(grid.get(s_c), 0),
             grobX(grid.get(e_c), 210), grobY(grid.get(e_c), 210),
             angle=90, shape= 1, square=FALSE,  curvature=0, open=TRUE, ncp=1,
             arrow=arr, gp=gpar(col=colors2[4], lwd=lwd * min_weight))
}

drawTransitions(2)
drawTransitions(1)



dev.off()

