
RunNo <- 33


# Enter Run number you want to graph

linetypes <- c(1,2,3)


par(mfrow=c(1,1), oma = c(1,1,1,1) + 0.2,  # CHANGE mfrow TO (2,2) for quad-plot 
   mar = c(3,3,1,1) + 1, mgp=c(1,1,0))
    
    
                                                                       
                                                       
matplot(cbind(X[RunNo,1:N], Y_old[RunNo,1:N],  Y[RunNo,1:N], b,0,1,b_trailing[RunNo,1:N], 
                Z_t_over_M_t[RunNo,1:N], Y_running_max[RunNo,1:N], X_double[RunNo, 1:N,s]), # tempstrike),  
        col=c("red", "green4",  "blue", "black", "black", "black",  "blue",  "darkorange", "darkorange", "black"), 
        type="l", lty=c(1, 1, 1, 2,2,2,1,2,1,1,1), lwd=0.1, xaxt = "n", yaxt = "n",
        main = "GBM reflected at trailing barrier (b x HWM)", xlab="Time steps", ylab="Price\n", 
        cex.main = 1.6, cex.axis = 1.6, cex.lab = 1.6, xlim=c((0/52*N),(52/52*N)) ,lim=c(0,1))  

text(4/52*N ,b-0.05,"fixed barrier", cex=1.4, col="black")
text(4/52*N ,b+0.09,"trailing barrier", cex=1.4, col="blue")
text(28/52*N ,2.7,"    GBM reflected at trailing barrier = U_t", cex=1.4, col= "blue")
text(27.5/52*N ,2.5,"    GBM reflected at fixed barrier = S_t", cex=1.4, col= "green4")
text(24.5/52*N ,2.3,"    Unrestricted GBM = X_t", cex=1.4, col= "red")
text(14/52*N ,0.23,"    Fractional drawdown, F_t = U_t / M_t    ", cex=1.4, col="darkorange")
text(20/52*N ,1.3,"    M_t", cex=1.4, col="darkorange")


axis(1, at=c(1,N), cex.axis = 1.4, labels = c("0", expression(italic("T"))))
axis(2, at=c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5), cex.axis = 1.4, labels = c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5), las=2)

