library(grid)
library(gridExtra)

# Layout Branches
pdf("matrices.pdf", width = 3.5, height=3)
grid.newpage()

names = c("Stick", "Branch", "Match", "Dark", "Merge")
heights = seq(.95,.1,-.8/5)
xs = rev(seq(0,.45,.45/5))
for(i in 1:length(names)) {
  pushViewport(viewport(x=xs[i], y = heights[i], width=0.5, height=.3, just=c("left","top"), name="vp1"))
  
  grid.rect(gp=gpar(fill="white", lwd=2))
  grid.grill()
  grid.text(names[i],
            x=unit(1, "npc")-unit(.5,"mm"),
            y=unit(1, "npc")-unit(1,"mm"), just=c("right", "top"),gp=gpar(fontstyle="bold", fontsize=9))
  popViewport()
}

dev.off()