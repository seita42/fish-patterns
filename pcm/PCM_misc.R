rateMat <- function(x) {
    rates <- matrix(c(    0, x$q12, x$q13, x$q14, x$q15, x$q16, x$q17, x$q18,
                      x$q21,     0, x$q23, x$q24, x$q25, x$q26, x$q27, x$q28,
                      x$q31, x$q32,     0, x$q34, x$q35, x$q36, x$q37, x$q38,
                      x$q41, x$q42, x$q43,     0, x$q45, x$q46, x$q47, x$q48,
                      x$q51, x$q52, x$q53, x$q54,     0, x$q56, x$q57, x$q58,
                      x$q61, x$q62, x$q63, x$q64, x$q65,     0, x$q67, x$q68,
                      x$q71, x$q72, x$q73, x$q74, x$q75, x$q76,     0, x$q78,
                      x$q81, x$q82, x$q83, x$q84, x$q85, x$q86, x$q87,     0),
                    nrow = 8, ncol = 8, byrow = TRUE)
    rownames(rates) <- colnames(rates) <- c("(0,0,0)", "(1,0,0)", "(0,1,0)", "(0,0,1)",
                                            "(1,1,0)", "(1,0,1)", "(0,1,1)", "(1,1,1)")
    return(rates)
}

rateMatMean <- function(x) {
    rates <- matrix(c(          0, mean(x$q12), mean(x$q13), mean(x$q14), mean(x$q15), mean(x$q16), mean(x$q17), mean(x$q18),
                      mean(x$q21),           0, mean(x$q23), mean(x$q24), mean(x$q25), mean(x$q26), mean(x$q27), mean(x$q28),
                      mean(x$q31), mean(x$q32),           0, mean(x$q34), mean(x$q35), mean(x$q36), mean(x$q37), mean(x$q38),
                      mean(x$q41), mean(x$q42), mean(x$q43),           0, mean(x$q45), mean(x$q46), mean(x$q47), mean(x$q48),
                      mean(x$q51), mean(x$q52), mean(x$q53), mean(x$q54),           0, mean(x$q56), mean(x$q57), mean(x$q58),
                      mean(x$q61), mean(x$q62), mean(x$q63), mean(x$q64), mean(x$q65),           0, mean(x$q67), mean(x$q68),
                      mean(x$q71), mean(x$q72), mean(x$q73), mean(x$q74), mean(x$q75), mean(x$q76),           0, mean(x$q78),
                      mean(x$q81), mean(x$q82), mean(x$q83), mean(x$q84), mean(x$q85), mean(x$q86), mean(x$q87),          0),
                    nrow = 8, ncol = 8, byrow = TRUE)
    rownames(rates) <- colnames(rates) <- c("(0,0,0)", "(1,0,0)", "(0,1,0)", "(0,0,1)",
                                            "(1,1,0)", "(1,0,1)", "(0,1,1)", "(1,1,1)")
    return(rates)
}

rateMatSd <- function(x) {
    rates <- matrix(c(          0, sd(x$q12), sd(x$q13), sd(x$q14), sd(x$q15), sd(x$q16), sd(x$q17), sd(x$q18),
                      sd(x$q21),           0, sd(x$q23), sd(x$q24), sd(x$q25), sd(x$q26), sd(x$q27), sd(x$q28),
                      sd(x$q31), sd(x$q32),           0, sd(x$q34), sd(x$q35), sd(x$q36), sd(x$q37), sd(x$q38),
                      sd(x$q41), sd(x$q42), sd(x$q43),           0, sd(x$q45), sd(x$q46), sd(x$q47), sd(x$q48),
                      sd(x$q51), sd(x$q52), sd(x$q53), sd(x$q54),           0, sd(x$q56), sd(x$q57), sd(x$q58),
                      sd(x$q61), sd(x$q62), sd(x$q63), sd(x$q64), sd(x$q65),           0, sd(x$q67), sd(x$q68),
                      sd(x$q71), sd(x$q72), sd(x$q73), sd(x$q74), sd(x$q75), sd(x$q76),           0, sd(x$q78),
                      sd(x$q81), sd(x$q82), sd(x$q83), sd(x$q84), sd(x$q85), sd(x$q86), sd(x$q87),          0),
                    nrow = 8, ncol = 8, byrow = TRUE)
    rownames(rates) <- colnames(rates) <- c("(0,0,0)", "(1,0,0)", "(0,1,0)", "(0,0,1)",
                                            "(1,1,0)", "(1,0,1)", "(0,1,1)", "(1,1,1)")
    return(rates)
}

rateMatMedian <- function(x) {
    rates <- matrix(c(          0, median(x$q12), median(x$q13), median(x$q14), median(x$q15), median(x$q16), median(x$q17), median(x$q18),
                      median(x$q21),           0, median(x$q23), median(x$q24), median(x$q25), median(x$q26), median(x$q27), median(x$q28),
                      median(x$q31), median(x$q32),           0, median(x$q34), median(x$q35), median(x$q36), median(x$q37), median(x$q38),
                      median(x$q41), median(x$q42), median(x$q43),           0, median(x$q45), median(x$q46), median(x$q47), median(x$q48),
                      median(x$q51), median(x$q52), median(x$q53), median(x$q54),           0, median(x$q56), median(x$q57), median(x$q58),
                      median(x$q61), median(x$q62), median(x$q63), median(x$q64), median(x$q65),           0, median(x$q67), median(x$q68),
                      median(x$q71), median(x$q72), median(x$q73), median(x$q74), median(x$q75), median(x$q76),           0, median(x$q78),
                      median(x$q81), median(x$q82), median(x$q83), median(x$q84), median(x$q85), median(x$q86), median(x$q87),          0),
                    nrow = 8, ncol = 8, byrow = TRUE)
    rownames(rates) <- colnames(rates) <- c("(0,0,0)", "(1,0,0)", "(0,1,0)", "(0,0,1)",
                                            "(1,1,0)", "(1,0,1)", "(0,1,1)", "(1,1,1)")
    return(rates)
}

ZrateMat <- function(x) {
    z <- x %>% summarize_at(vars(starts_with("q")), list(~ sum(.==0.0)/n()))
    Zrates <- matrix(c(    0, z$q12, z$q13, z$q14, z$q15, z$q16, z$q17, z$q18,
                       z$q21,     0, z$q23, z$q24, z$q25, z$q26, z$q27, z$q28,
                       z$q31, z$q32,     0, z$q34, z$q35, z$q36, z$q37, z$q38,
                       z$q41, z$q42, z$q43,     0, z$q45, z$q46, z$q47, z$q48,
                       z$q51, z$q52, z$q53, z$q54,     0, z$q56, z$q57, z$q58,
                       z$q61, z$q62, z$q63, z$q64, z$q65,     0, z$q67, z$q68,
                       z$q71, z$q72, z$q73, z$q74, z$q75, z$q76,     0, z$q78,
                       z$q81, z$q82, z$q83, z$q84, z$q85, z$q86, z$q87,     0),
                     nrow = 8, ncol = 8, byrow = TRUE)
    rownames(Zrates) <- colnames(Zrates) <- c("(0,0,0)", "(1,0,0)", "(0,1,0)", "(0,0,1)",
                                              "(1,1,0)", "(1,0,1)", "(0,1,1)", "(1,1,1)")
    return(Zrates)
}

plotRates <- function(x, ...) {
    x[is.na(x)] <- 0
    # print(x)

    if (hasArg(signif)) 
        signif <- list(...)$signif
    else signif <- 3
    if (hasArg(main)) 
        main <- list(...)$main
    else main <- NULL
    if (hasArg(cex.main)) 
        cex.main <- list(...)$cex.main
    else cex.main <- 1.2
    if (hasArg(sub)) 
        sub <- list(...)$sub
    else sub <- NULL
    if (hasArg(cex.sub)) 
        cex.sub <- list(...)$cex.sub
    else cex.sub <- 1
    if (hasArg(cex.traits)) 
        cex.traits <- list(...)$cex.traits
    else cex.traits <- 0.9
    if (hasArg(cex.rates)) 
        cex.rates <- list(...)$cex.rates
    else cex.rates <- 0.8
    if (hasArg(lwd.by.rate))
        lwd.by.rate <- list(...)$lwd.by.rate
    else lwd.by.rate <- FALSE
    if (lwd.by.rate) {
        rates <- c(x[x > 0])
        LWD <- round(x/min(rates))
        if (hasArg(max.lwd)) 
            max.lwd <- list(...)$max.lwd
        else max.lwd <- 10
        LWD[LWD > max.lwd] <- max.lwd
    }
    else LWD <- matrix(2, nrow(x), ncol(x))
    
    plot.new()
    par(mar = c(1.1, 2.1, 3.1, 2.1))
    plot.window(xlim = c(0, 2), ylim = c(0, 2), asp = 1)
    mtext(
        if (!is.null(main))
            main
        else
            "Rates"
        , side = 3, adj = 0, line = 1.2,
        cex = cex.main)
    text(x = 0.0, y = 2.0, adj = 0,
        if (!is.null(sub))
            sub
        else
            "(X, Y, Z):"
        , cex = cex.sub)
    
    px <- c(0.2, 0.65, 1.35, 1.8)
    py <- c(1.8, 1.35, 0.65, 0.2)
    adx <- 0.15
    ady <- 0.1
    adx2 <- 0.1
    ady2 <- 0.1
    amx <- 0.05
    amy <- 0.05
    aix <- 0.025
    aiy <- 0.025
    tdx <- 0.08
    tdy <- 0.08
    
    px1 <- px[1]
    px2 <- px[1]
    px3 <- px[2]
    px4 <- px[4]
    px5 <- px[2]
    px6 <- px[4]
    px7 <- px[3]
    px8 <- px[3]
    
    py1 <- py[2]
    py2 <- py[3]
    py3 <- py[1]
    py4 <- py[2]
    py5 <- py[4]
    py6 <- py[3]
    py7 <- py[1]
    py8 <- py[4]
        

    # Vertical
    
    # 000(1) <-> 100(2)
    arrows(x0 = px1-aix, y0 = py2+ady, x1 = px1-aix, y1 = py1-ady,
           lwd = max(LWD[2, 1], 1),
           lty = if (LWD[2, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px1+aix, y0 = py1-ady, x1 = px1+aix, y1 = py2+ady,
           lwd = max(LWD[1, 2], 1),
           lty = if (LWD[1, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 001(4) <-> 101(6)
    arrows(x0 = px4-aix, y0 = py6+ady, x1 = px4-aix, y1 = py4-ady,
           lwd = max(LWD[6, 4], 1),
           lty = if (LWD[6, 4] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px4+aix, y0 = py4-ady, x1 = px4+aix, y1 = py6+ady,
           lwd = max(LWD[4, 6], 1),
           lty = if (LWD[4, 6] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 010(3) <-> 110(5)
    arrows(x0 = px3-aix, y0 = py5+ady, x1 = px3-aix, y1 = py3-ady,
           lwd = max(LWD[5, 3], 1),
           lty = if (LWD[5, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px3+aix, y0 = py3-ady, x1 = px3+aix, y1 = py5+ady,
           lwd = max(LWD[3, 5], 1),
           lty = if (LWD[3, 5] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 011(7) <-> 111(8)
    arrows(x0 = px7-aix, y0 = py8+ady, x1 = px7-aix, y1 = py7-ady,
           lwd = max(LWD[8, 7], 1),
           lty = if (LWD[8, 7] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px7+aix, y0 = py7-ady, x1 = px7+aix, y1 = py8+ady,
           lwd = max(LWD[7, 8], 1),
           lty = if (LWD[7, 8] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)


    # Horizontal
    
    # 000(1) <-> 001(4)
    arrows(x0 = px1+adx, y0 = py1+aiy, x1 = px4-adx, y1 = py1+aiy,
           lwd = max(LWD[1, 4], 1),
           lty = if (LWD[1, 4] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px4-adx, y0 = py1-aiy, x1 = px1+adx, y1 = py1-aiy,
           lwd = max(LWD[4, 1], 1),
           lty = if (LWD[4, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 100(2) <-> 101(6)
    arrows(x0 = px2+adx, y0 = py2+aiy, x1 = px6-adx, y1 = py2+aiy,
           lwd = max(LWD[2, 6], 1),
           lty = if (LWD[2, 6] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px6-adx, y0 = py2-aiy, x1 = px2+adx, y1 = py2-aiy,
           lwd = max(LWD[6, 2], 1),
           lty = if (LWD[6, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 010(3) <-> 011(7)
    arrows(x0 = px3+adx, y0 = py3+aiy, x1 = px7-adx, y1 = py3+aiy,
           lwd = max(LWD[3, 7], 1),
           lty = if (LWD[3, 7] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px7-adx, y0 = py3-aiy, x1 = px3+adx, y1 = py3-aiy,
           lwd = max(LWD[7, 3], 1),
           lty = if (LWD[7, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 110(5) <-> 111(8)
    arrows(x0 = px5+adx, y0 = py5+aiy, x1 = px8-adx, y1 = py5+aiy,
           lwd = max(LWD[5, 8], 1),
           lty = if (LWD[5, 8] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px8-adx, y0 = py5-aiy, x1 = px5+adx, y1 = py5-aiy,
           lwd = max(LWD[8, 5], 1),
           lty = if (LWD[8, 5] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)


    # Diagonal
    
    # 000(1) <-> 010(3)
    arrows(x0 = px1+(adx2-aiy)*0.7-amx, y0 = py1+(ady2+aix)*0.7+amy,
           x1 = px3-(adx2+aiy)*0.7-amx, y1 = py3-(ady2-aix)*0.7+amy,
           lwd = max(LWD[1, 3], 1),
           lty = if (LWD[1, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px3-(adx2-aiy)*0.7-amx, y0 = py3-(ady2+aix)*0.7+amy,
           x1 = px1+(adx2+aiy)*0.7-amx, y1 = py1+(ady2-aix)*0.7+amy,
           lwd = max(LWD[3, 1], 1),
           lty = if (LWD[3, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 001(4) <-> 011(7)
    arrows(x0 = px4-(adx2+aiy)*0.7+amx, y0 = py4+(ady2-aix)*0.7+amy,
           x1 = px7+(adx2-aiy)*0.7+amx, y1 = py7-(ady2+aix)*0.7+amy,
           lwd = max(LWD[4, 7], 1),
           lty = if (LWD[4, 7] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px7+(adx2+aiy)*0.7+amx, y0 = py7-(ady2-aix)*0.7+amy,
           x1 = px4-(adx2-aiy)*0.7+amx, y1 = py4+(ady2+aix)*0.7+amy,
           lwd = max(LWD[7, 4], 1),
           lty = if (LWD[7, 4] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 100(2) <-> 110(5)
    arrows(x0 = px2+(adx2+aiy)*0.7-amx, y0 = py2-(ady2-aix)*0.7-amy,
           x1 = px5-(adx2-aiy)*0.7-amx, y1 = py5+(ady2+aix)*0.7-amy,
           lwd = max(LWD[2, 5], 1),
           lty = if (LWD[2, 5] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px5-(adx2+aiy)*0.7-amx, y0 = py5+(ady2-aix)*0.7-amy,
           x1 = px2+(adx2-aiy)*0.7-amx, y1 = py2-(ady2+aix)*0.7-amy,
           lwd = max(LWD[5, 2], 1),
           lty = if (LWD[5, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 101(6) <-> 111(8)
    arrows(x0 = px6-(adx2-aiy)*0.7+amx, y0 = py6-(ady2+aix)*0.7-amy,
           x1 = px8+(adx2+aiy)*0.7+amx, y1 = py8+(ady2-aix)*0.7-amy,
           lwd = max(LWD[6, 8], 1),
           lty = if (LWD[6, 8] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px8+(adx2-aiy)*0.7+amx, y0 = py8+(ady2+aix)*0.7-amy,
           x1 = px6-(adx2+aiy)*0.7+amx, y1 = py6-(ady2-aix)*0.7-amy,
           lwd = max(LWD[8, 6], 1),
           lty = if (LWD[8, 6] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    text(x = px1, y = py1, rownames(x)[1], cex = cex.traits)
    text(x = px2, y = py2, rownames(x)[2], cex = cex.traits)
    text(x = px3, y = py3, rownames(x)[3], cex = cex.traits)
    text(x = px4, y = py4, rownames(x)[4], cex = cex.traits)
    text(x = px5, y = py5, rownames(x)[5], cex = cex.traits)
    text(x = px6, y = py6, rownames(x)[6], cex = cex.traits)
    text(x = px7, y = py7, rownames(x)[7], cex = cex.traits)
    text(x = px8, y = py8, rownames(x)[8], cex = cex.traits)

    text(x = (px3+px7)/2.0, y = py3+tdy, round(x[3, 7], signif), cex = cex.rates)
    text(x = (px3+px7)/2.0, y = py3-tdy, round(x[7, 3], signif), cex = cex.rates)
    text(x = (px1+px4)/2.0, y = py1+tdy, round(x[1, 4], signif), cex = cex.rates)
    text(x = (px1+px4)/2.0, y = py1-tdy, round(x[4, 1], signif), cex = cex.rates)
    text(x = (px2+px6)/2.0, y = py2+tdy, round(x[2, 6], signif), cex = cex.rates)
    text(x = (px2+px6)/2.0, y = py2-tdy, round(x[6, 2], signif), cex = cex.rates)
    text(x = (px5+px8)/2.0, y = py5+tdy, round(x[5, 8], signif), cex = cex.rates)
    text(x = (px5+px8)/2.0, y = py5-tdy, round(x[8, 5], signif), cex = cex.rates)

    text(x = px1-tdx, y = (py1+py2)/2.0, round(x[2, 1], signif), cex = cex.rates, srt = 90)
    text(x = px1+tdx, y = (py1+py2)/2.0, round(x[1, 2], signif), cex = cex.rates, srt = 90)
    text(x = px3-tdx, y = (py3+py5)/2.0, round(x[5, 3], signif), cex = cex.rates, srt = 90)
    text(x = px3+tdx, y = (py3+py5)/2.0, round(x[3, 5], signif), cex = cex.rates, srt = 90)
    text(x = px7-tdx, y = (py7+py8)/2.0, round(x[8, 7], signif), cex = cex.rates, srt = 90)
    text(x = px7+tdx, y = (py7+py8)/2.0, round(x[7, 8], signif), cex = cex.rates, srt = 90)
    text(x = px4-tdx, y = (py4+py6)/2.0, round(x[6, 4], signif), cex = cex.rates, srt = 90)
    text(x = px4+tdx, y = (py4+py6)/2.0, round(x[4, 6], signif), cex = cex.rates, srt = 90)

    text(x = (px1+px3)/2.0-tdx*0.7-amx, y = (py1+py3)/2.0+tdy*0.7+amy, round(x[1, 3], signif), cex = cex.rates, srt = 45)
    text(x = (px1+px3)/2.0+tdx*0.7-amx, y = (py1+py3)/2.0-tdy*0.7+amy, round(x[3, 1], signif), cex = cex.rates, srt = 45)
    text(x = (px4+px7)/2.0-tdx*0.7+amx, y = (py4+py7)/2.0-tdy*0.7+amy, round(x[4, 7], signif), cex = cex.rates, srt = -45)
    text(x = (px4+px7)/2.0+tdx*0.7+amx, y = (py4+py7)/2.0+tdy*0.7+amy, round(x[7, 4], signif), cex = cex.rates, srt = -45)
    text(x = (px2+px5)/2.0+tdx*0.7-amx, y = (py2+py5)/2.0+tdy*0.7-amy, round(x[2, 5], signif), cex = cex.rates, srt = -45)
    text(x = (px2+px5)/2.0-tdx*0.7-amx, y = (py2+py5)/2.0-tdy*0.7-amy, round(x[5, 2], signif), cex = cex.rates, srt = -45)
    text(x = (px6+px8)/2.0+tdx*0.7+amx, y = (py6+py8)/2.0-tdy*0.7-amy, round(x[6, 8], signif), cex = cex.rates, srt = 45)
    text(x = (px6+px8)/2.0-tdx*0.7+amx, y = (py6+py8)/2.0+tdy*0.7-amy, round(x[8, 6], signif), cex = cex.rates, srt = 45)

}

plotRatesInd <- function(x, ...) {
    x[is.na(x)] <- 0
    # print(x)

    if (hasArg(signif)) 
        signif <- list(...)$signif
    else signif <- 3
    if (hasArg(main)) 
        main <- list(...)$main
    else main <- NULL
    if (hasArg(cex.main)) 
        cex.main <- list(...)$cex.main
    else cex.main <- 1.2
    if (hasArg(sub)) 
        sub <- list(...)$sub
    else sub <- NULL
    if (hasArg(cex.sub)) 
        cex.sub <- list(...)$cex.sub
    else cex.sub <- 1
    if (hasArg(cex.traits)) 
        cex.traits <- list(...)$cex.traits
    else cex.traits <- 0.9
    if (hasArg(cex.rates)) 
        cex.rates <- list(...)$cex.rates
    else cex.rates <- 0.8
    if (hasArg(lwd.by.rate))
        lwd.by.rate <- list(...)$lwd.by.rate
    else lwd.by.rate <- FALSE
    if (lwd.by.rate) {
        rates <- c(x[x > 0])
        LWD <- round(x/min(rates))
        if (hasArg(max.lwd)) 
            max.lwd <- list(...)$max.lwd
        else max.lwd <- 10
        LWD[LWD > max.lwd] <- max.lwd
    }
    else LWD <- matrix(2, nrow(x), ncol(x))
    
    options(repr.plot.width = 8, repr.plot.height = 4)
    plot.new()
    par(mar = c(1.1, 2.1, 3.1, 2.1))
    plot.window(xlim = c(0, 2), ylim = c(0, 1), asp = 1)
    mtext(
        if (!is.null(main))
            main
        else
            "Rates"
        , side = 3, adj = 0, line = 1.2,
        cex = cex.main)
    text(x = 0.0, y = 1.0, adj = 0,
        if (!is.null(sub))
            sub
        else
            "(X, Y, Z):"
        , cex = cex.sub)
    
    px <- c(0.4, 1.0, 1.6)
    py <- c(0.8, 0.2)
    adx <- 0.15
    ady <- 0.1
    adx2 <- 0.1
    ady2 <- 0.1
    amx <- 0.05
    amy <- 0.05
    aix <- 0.025
    aiy <- 0.025
    tdx <- 0.08
    tdy <- 0.08
    
    px1 <- px[1]
    px2 <- px[1]
    px3 <- px[2]
    px4 <- px[2]
    px5 <- px[3]
    px6 <- px[3]
    
    py1 <- py[1]
    py2 <- py[2]
    py3 <- py[1]
    py4 <- py[2]
    py5 <- py[1]
    py6 <- py[2]
        

    # Vertical
    
    # (0,*,*)(1347) <-> (1,*,*)(2568)
    arrows(x0 = px2-aix, y0 = py2+ady, x1 = px1-aix, y1 = py1-ady,
           lwd = max(LWD[2, 1], 1),
           lty = if (LWD[2, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px1+aix, y0 = py1-ady, x1 = px2+aix, y1 = py2+ady,
           lwd = max(LWD[1, 2], 1),
           lty = if (LWD[1, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # (*,0,*)(1246) <-> (*,1,*)(3578)
    arrows(x0 = px4-aix, y0 = py4+ady, x1 = px3-aix, y1 = py3-ady,
           lwd = max(LWD[3, 1], 1),
           lty = if (LWD[3, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px3+aix, y0 = py3-ady, x1 = px4+aix, y1 = py4+ady,
           lwd = max(LWD[1, 3], 1),
           lty = if (LWD[1, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # (*,*,0)(1235) <-> (*,*,1)(4678)
    arrows(x0 = px6-aix, y0 = py6+ady, x1 = px5-aix, y1 = py5-ady,
           lwd = max(LWD[4, 1], 1),
           lty = if (LWD[4, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px5+aix, y0 = py5-ady, x1 = px6+aix, y1 = py6+ady,
           lwd = max(LWD[1, 4], 1),
           lty = if (LWD[1, 4] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    
    ratetext <- function(x) {
        return(round(x, signif))
    }

    text(x = px1, y = py1, "(0,*,*)", cex = cex.traits)
    text(x = px2, y = py2, "(1,*,*)", cex = cex.traits)
    text(x = px3, y = py3, "(*,0,*)", cex = cex.traits)
    text(x = px4, y = py4, "(*,1,*)", cex = cex.traits)
    text(x = px5, y = py5, "(*,*,0)", cex = cex.traits)
    text(x = px6, y = py6, "(*,*,1)", cex = cex.traits)

    text(x = px1-tdx, y = (py1+py2)/2.0, ratetext(x[2, 1]), cex = cex.rates, srt = 90)
    text(x = px1+tdx, y = (py1+py2)/2.0, ratetext(x[1, 2]), cex = cex.rates, srt = 90)
    text(x = px3-tdx, y = (py3+py4)/2.0, ratetext(x[3, 1]), cex = cex.rates, srt = 90)
    text(x = px3+tdx, y = (py3+py4)/2.0, ratetext(x[1, 3]), cex = cex.rates, srt = 90)
    text(x = px5-tdx, y = (py5+py6)/2.0, ratetext(x[4, 1]), cex = cex.rates, srt = 90)
    text(x = px5+tdx, y = (py5+py6)/2.0, ratetext(x[1, 4]), cex = cex.rates, srt = 90)

}

plotRatesIndZ <- function(x, ...) {
    x[is.na(x)] <- 0
    # print(x)

    if (hasArg(signif)) 
        signif <- list(...)$signif
    else signif <- 3
    if (hasArg(main)) 
        main <- list(...)$main
    else main <- NULL
    if (hasArg(cex.main)) 
        cex.main <- list(...)$cex.main
    else cex.main <- 1.2
    if (hasArg(sub)) 
        sub <- list(...)$sub
    else sub <- NULL
    if (hasArg(cex.sub)) 
        cex.sub <- list(...)$cex.sub
    else cex.sub <- 1
    if (hasArg(cex.traits)) 
        cex.traits <- list(...)$cex.traits
    else cex.traits <- 0.9
    if (hasArg(cex.rates)) 
        cex.rates <- list(...)$cex.rates
    else cex.rates <- 0.8
    if (hasArg(lwd.by.rate))
        lwd.by.rate <- list(...)$lwd.by.rate
    else lwd.by.rate <- FALSE
    if (lwd.by.rate) {
        rates <- c(x[x > 0])
        LWD <- round(x/min(rates))
        if (hasArg(max.lwd)) 
            max.lwd <- list(...)$max.lwd
        else max.lwd <- 10
        LWD[LWD > max.lwd] <- max.lwd
    }
    else LWD <- matrix(2, nrow(x), ncol(x))
    
    options(repr.plot.width = 8, repr.plot.height = 4)
    plot.new()
    par(mar = c(1.1, 2.1, 3.1, 2.1))
    plot.window(xlim = c(0, 2), ylim = c(0, 1), asp = 1)
    mtext(
        if (!is.null(main))
            main
        else
            "Rates"
        , side = 3, adj = 0, line = 1.2,
        cex = cex.main)
    text(x = 0.0, y = 1.0, adj = 0,
        if (!is.null(sub))
            sub
        else
            "(X, Y, Z):"
        , cex = cex.sub)
    
    px <- c(0.2, 1.2, 1.8)
    py <- c(0.8, 0.2)
    adx <- 0.15
    ady <- 0.1
    adx2 <- 0.1
    ady2 <- 0.1
    amx <- 0.05
    amy <- 0.05
    aix <- 0.025
    aiy <- 0.025
    tdx <- 0.08
    tdy <- 0.08
    
    px1 <- px[1]  # (0,0,*)
    px2 <- px[1]  # (1,0,*)
    px3 <- px[2]  # (0,1,*)
    px4 <- px[2]  # (1,1,*)
    px5 <- px[3]  # (*,*,0)
    px6 <- px[3]  # (*,*,1)
    
    py1 <- py[1]
    py2 <- py[2]
    py3 <- py[1]
    py4 <- py[2]
    py5 <- py[1]
    py6 <- py[2]
        

    # Vertical
    
    # (0,0,*)(14) <-> (1,0,*)(26)
    arrows(x0 = px2-aix, y0 = py2+ady, x1 = px1-aix, y1 = py1-ady,
           lwd = max(LWD[2, 1], 1),
           lty = if (LWD[2, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px1+aix, y0 = py1-ady, x1 = px2+aix, y1 = py2+ady,
           lwd = max(LWD[1, 2], 1),
           lty = if (LWD[1, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # (0,1,*)(37) <-> (1,1,*)(58)
    arrows(x0 = px4-aix, y0 = py4+ady, x1 = px3-aix, y1 = py3-ady,
           lwd = max(LWD[5, 3], 1),
           lty = if (LWD[5, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px3+aix, y0 = py3-ady, x1 = px4+aix, y1 = py4+ady,
           lwd = max(LWD[3, 5], 1),
           lty = if (LWD[3, 5] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # (*,*,0)(1235) <-> (*,*,1)(4678)
    arrows(x0 = px6-aix, y0 = py6+ady, x1 = px5-aix, y1 = py5-ady,
           lwd = max(LWD[4, 1], 1),
           lty = if (LWD[4, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px5+aix, y0 = py5-ady, x1 = px6+aix, y1 = py6+ady,
           lwd = max(LWD[1, 4], 1),
           lty = if (LWD[1, 4] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # Horizontal
    
    # (0,0,*)(14) <-> (0,1,*)(37)
    arrows(x0 = px1+adx, y0 = py1+aiy, x1 = px3-adx, y1 = py3+aiy,
           lwd = max(LWD[1, 3], 1),
           lty = if (LWD[1, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px3-adx, y0 = py3-aiy, x1 = px1+adx, y1 = py1-aiy,
           lwd = max(LWD[3, 1], 1),
           lty = if (LWD[3, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # (1,0,*)(26) <-> (1,1,*)(58)
    arrows(x0 = px2+adx, y0 = py2+aiy, x1 = px4-adx, y1 = py4+aiy,
           lwd = max(LWD[2, 5], 1),
           lty = if (LWD[2, 5] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px4-adx, y0 = py4-aiy, x1 = px2+adx, y1 = py2-aiy,
           lwd = max(LWD[5, 2], 1),
           lty = if (LWD[5, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    
    ratetext <- function(x) {
        return(round(x, signif))
    }

    text(x = px1, y = py1, "(0,0,*)", cex = cex.traits)
    text(x = px2, y = py2, "(1,0,*)", cex = cex.traits)
    text(x = px3, y = py3, "(0,1,*)", cex = cex.traits)
    text(x = px4, y = py4, "(1,1,*)", cex = cex.traits)
    text(x = px5, y = py5, "(*,*,0)", cex = cex.traits)
    text(x = px6, y = py6, "(*,*,1)", cex = cex.traits)

    text(x = px1-tdx, y = (py1+py2)/2.0, ratetext(x[2, 1]), cex = cex.rates, srt = 90)
    text(x = px1+tdx, y = (py1+py2)/2.0, ratetext(x[1, 2]), cex = cex.rates, srt = 90)
    text(x = px3-tdx, y = (py3+py4)/2.0, ratetext(x[5, 3]), cex = cex.rates, srt = 90)
    text(x = px3+tdx, y = (py3+py4)/2.0, ratetext(x[3, 5]), cex = cex.rates, srt = 90)
    text(x = px5-tdx, y = (py5+py6)/2.0, ratetext(x[4, 1]), cex = cex.rates, srt = 90)
    text(x = px5+tdx, y = (py5+py6)/2.0, ratetext(x[1, 4]), cex = cex.rates, srt = 90)
    
    text(x = (px1+px3)/2.0, y = py1+tdy, ratetext(x[1, 3]), cex = cex.rates)
    text(x = (px1+px3)/2.0, y = py1-tdy, ratetext(x[3, 1]), cex = cex.rates)
    text(x = (px2+px4)/2.0, y = py2+tdy, ratetext(x[2, 5]), cex = cex.rates)
    text(x = (px2+px4)/2.0, y = py2-tdy, ratetext(x[5, 2]), cex = cex.rates)

}

plotRatesZ <- function(x, z, ...) {
    x[is.na(x)] <- 0
    # print(x)
    z[is.na(z)] <- 0
    # print(z)

    if (hasArg(signif)) 
        signif <- list(...)$signif
    else signif <- 3
    if (hasArg(main)) 
        main <- list(...)$main
    else main <- NULL
    if (hasArg(cex.main)) 
        cex.main <- list(...)$cex.main
    else cex.main <- 1.2
    if (hasArg(sub)) 
        sub <- list(...)$sub
    else sub <- NULL
    if (hasArg(cex.sub)) 
        cex.sub <- list(...)$cex.sub
    else cex.sub <- 1
    if (hasArg(cex.traits)) 
        cex.traits <- list(...)$cex.traits
    else cex.traits <- 0.9
    if (hasArg(cex.rates)) 
        cex.rates <- list(...)$cex.rates
    else cex.rates <- 0.8
    if (hasArg(lwd.by.rate))
        lwd.by.rate <- list(...)$lwd.by.rate
    else lwd.by.rate <- FALSE
    if (lwd.by.rate) {
        rates <- c(x[x > 0])
        LWD <- round(x/min(rates))
        if (hasArg(max.lwd)) 
            max.lwd <- list(...)$max.lwd
        else max.lwd <- 10
        LWD[LWD > max.lwd] <- max.lwd
    }
    else LWD <- matrix(2, nrow(x), ncol(x))
    
    options(repr.plot.width = 8, repr.plot.height = 8)
    plot.new()
    par(mar = c(1.1, 2.1, 3.1, 2.1))
    plot.window(xlim = c(0, 2), ylim = c(0, 2), asp = 1)
    mtext(
        if (!is.null(main))
            main
        else
            "Rates"
        , side = 3, adj = 0, line = 1.2,
        cex = cex.main)
    text(x = 0.0, y = 2.0, adj = 0,
        if (!is.null(sub))
            sub
        else
            "(X, Y, Z):"
        , cex = cex.sub)
    
    px <- c(0.2, 0.65, 1.35, 1.8)
    py <- c(1.8, 1.35, 0.65, 0.2)
    adx <- 0.15
    ady <- 0.1
    adx2 <- 0.1
    ady2 <- 0.1
    amx <- 0.05
    amy <- 0.05
    aix <- 0.025
    aiy <- 0.025
    tdx <- 0.08
    tdy <- 0.08
    
    px1 <- px[1]
    px2 <- px[1]
    px3 <- px[2]
    px4 <- px[4]
    px5 <- px[2]
    px6 <- px[4]
    px7 <- px[3]
    px8 <- px[3]
    
    py1 <- py[2]
    py2 <- py[3]
    py3 <- py[1]
    py4 <- py[2]
    py5 <- py[4]
    py6 <- py[3]
    py7 <- py[1]
    py8 <- py[4]
        

    # Vertical
    
    # 000(1) <-> 100(2)
    arrows(x0 = px1-aix, y0 = py2+ady, x1 = px1-aix, y1 = py1-ady,
           lwd = max(LWD[2, 1], 1),
           lty = if (LWD[2, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px1+aix, y0 = py1-ady, x1 = px1+aix, y1 = py2+ady,
           lwd = max(LWD[1, 2], 1),
           lty = if (LWD[1, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 001(4) <-> 101(6)
    arrows(x0 = px4-aix, y0 = py6+ady, x1 = px4-aix, y1 = py4-ady,
           lwd = max(LWD[6, 4], 1),
           lty = if (LWD[6, 4] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px4+aix, y0 = py4-ady, x1 = px4+aix, y1 = py6+ady,
           lwd = max(LWD[4, 6], 1),
           lty = if (LWD[4, 6] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 010(3) <-> 110(5)
    arrows(x0 = px3-aix, y0 = py5+ady, x1 = px3-aix, y1 = py3-ady,
           lwd = max(LWD[5, 3], 1),
           lty = if (LWD[5, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px3+aix, y0 = py3-ady, x1 = px3+aix, y1 = py5+ady,
           lwd = max(LWD[3, 5], 1),
           lty = if (LWD[3, 5] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 011(7) <-> 111(8)
    arrows(x0 = px7-aix, y0 = py8+ady, x1 = px7-aix, y1 = py7-ady,
           lwd = max(LWD[8, 7], 1),
           lty = if (LWD[8, 7] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px7+aix, y0 = py7-ady, x1 = px7+aix, y1 = py8+ady,
           lwd = max(LWD[7, 8], 1),
           lty = if (LWD[7, 8] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)


    # Horizontal
    
    # 000(1) <-> 001(4)
    arrows(x0 = px1+adx, y0 = py1+aiy, x1 = px4-adx, y1 = py1+aiy,
           lwd = max(LWD[1, 4], 1),
           lty = if (LWD[1, 4] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px4-adx, y0 = py1-aiy, x1 = px1+adx, y1 = py1-aiy,
           lwd = max(LWD[4, 1], 1),
           lty = if (LWD[4, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 100(2) <-> 101(6)
    arrows(x0 = px2+adx, y0 = py2+aiy, x1 = px6-adx, y1 = py2+aiy,
           lwd = max(LWD[2, 6], 1),
           lty = if (LWD[2, 6] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px6-adx, y0 = py2-aiy, x1 = px2+adx, y1 = py2-aiy,
           lwd = max(LWD[6, 2], 1),
           lty = if (LWD[6, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 010(3) <-> 011(7)
    arrows(x0 = px3+adx, y0 = py3+aiy, x1 = px7-adx, y1 = py3+aiy,
           lwd = max(LWD[3, 7], 1),
           lty = if (LWD[3, 7] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px7-adx, y0 = py3-aiy, x1 = px3+adx, y1 = py3-aiy,
           lwd = max(LWD[7, 3], 1),
           lty = if (LWD[7, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 110(5) <-> 111(8)
    arrows(x0 = px5+adx, y0 = py5+aiy, x1 = px8-adx, y1 = py5+aiy,
           lwd = max(LWD[5, 8], 1),
           lty = if (LWD[5, 8] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px8-adx, y0 = py5-aiy, x1 = px5+adx, y1 = py5-aiy,
           lwd = max(LWD[8, 5], 1),
           lty = if (LWD[8, 5] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)


    # Diagonal
    
    # 000(1) <-> 010(3)
    arrows(x0 = px1+(adx2-aiy)*0.7-amx, y0 = py1+(ady2+aix)*0.7+amy,
           x1 = px3-(adx2+aiy)*0.7-amx, y1 = py3-(ady2-aix)*0.7+amy,
           lwd = max(LWD[1, 3], 1),
           lty = if (LWD[1, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px3-(adx2-aiy)*0.7-amx, y0 = py3-(ady2+aix)*0.7+amy,
           x1 = px1+(adx2+aiy)*0.7-amx, y1 = py1+(ady2-aix)*0.7+amy,
           lwd = max(LWD[3, 1], 1),
           lty = if (LWD[3, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 001(4) <-> 011(7)
    arrows(x0 = px4-(adx2+aiy)*0.7+amx, y0 = py4+(ady2-aix)*0.7+amy,
           x1 = px7+(adx2-aiy)*0.7+amx, y1 = py7-(ady2+aix)*0.7+amy,
           lwd = max(LWD[4, 7], 1),
           lty = if (LWD[4, 7] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px7+(adx2+aiy)*0.7+amx, y0 = py7-(ady2-aix)*0.7+amy,
           x1 = px4-(adx2-aiy)*0.7+amx, y1 = py4+(ady2+aix)*0.7+amy,
           lwd = max(LWD[7, 4], 1),
           lty = if (LWD[7, 4] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 100(2) <-> 110(5)
    arrows(x0 = px2+(adx2+aiy)*0.7-amx, y0 = py2-(ady2-aix)*0.7-amy,
           x1 = px5-(adx2-aiy)*0.7-amx, y1 = py5+(ady2+aix)*0.7-amy,
           lwd = max(LWD[2, 5], 1),
           lty = if (LWD[2, 5] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px5-(adx2+aiy)*0.7-amx, y0 = py5+(ady2-aix)*0.7-amy,
           x1 = px2+(adx2-aiy)*0.7-amx, y1 = py2-(ady2+aix)*0.7-amy,
           lwd = max(LWD[5, 2], 1),
           lty = if (LWD[5, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # 101(6) <-> 111(8)
    arrows(x0 = px6-(adx2-aiy)*0.7+amx, y0 = py6-(ady2+aix)*0.7-amy,
           x1 = px8+(adx2+aiy)*0.7+amx, y1 = py8+(ady2-aix)*0.7-amy,
           lwd = max(LWD[6, 8], 1),
           lty = if (LWD[6, 8] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px8+(adx2-aiy)*0.7+amx, y0 = py8+(ady2+aix)*0.7-amy,
           x1 = px6-(adx2+aiy)*0.7+amx, y1 = py6-(ady2-aix)*0.7-amy,
           lwd = max(LWD[8, 6], 1),
           lty = if (LWD[8, 6] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    
    ratetext <- function(x, z) {
        return(paste(round(x, signif), " (z = ", round(z, signif), ")", sep=""))
    }

    text(x = px1, y = py1, rownames(x)[1], cex = cex.traits)
    text(x = px2, y = py2, rownames(x)[2], cex = cex.traits)
    text(x = px3, y = py3, rownames(x)[3], cex = cex.traits)
    text(x = px4, y = py4, rownames(x)[4], cex = cex.traits)
    text(x = px5, y = py5, rownames(x)[5], cex = cex.traits)
    text(x = px6, y = py6, rownames(x)[6], cex = cex.traits)
    text(x = px7, y = py7, rownames(x)[7], cex = cex.traits)
    text(x = px8, y = py8, rownames(x)[8], cex = cex.traits)

    text(x = (px3+px7)/2.0, y = py3+tdy, ratetext(x[3, 7], z[3, 7]), cex = cex.rates)
    text(x = (px3+px7)/2.0, y = py3-tdy, ratetext(x[7, 3], z[7, 3]), cex = cex.rates)
    text(x = (px1+px4)/2.0, y = py1+tdy, ratetext(x[1, 4], z[1, 4]), cex = cex.rates)
    text(x = (px1+px4)/2.0, y = py1-tdy, ratetext(x[4, 1], z[4, 1]), cex = cex.rates)
    text(x = (px2+px6)/2.0, y = py2+tdy, ratetext(x[2, 6], z[2, 6]), cex = cex.rates)
    text(x = (px2+px6)/2.0, y = py2-tdy, ratetext(x[6, 2], z[6, 2]), cex = cex.rates)
    text(x = (px5+px8)/2.0, y = py5+tdy, ratetext(x[5, 8], z[5, 8]), cex = cex.rates)
    text(x = (px5+px8)/2.0, y = py5-tdy, ratetext(x[8, 5], z[8, 5]), cex = cex.rates)

    text(x = px1-tdx, y = (py1+py2)/2.0, ratetext(x[2, 1], z[2, 1]), cex = cex.rates, srt = 90)
    text(x = px1+tdx, y = (py1+py2)/2.0, ratetext(x[1, 2], z[1, 2]), cex = cex.rates, srt = 90)
    text(x = px3-tdx, y = (py3+py5)/2.0, ratetext(x[5, 3], z[5, 3]), cex = cex.rates, srt = 90)
    text(x = px3+tdx, y = (py3+py5)/2.0, ratetext(x[3, 5], z[3, 5]), cex = cex.rates, srt = 90)
    text(x = px7-tdx, y = (py7+py8)/2.0, ratetext(x[8, 7], z[8, 7]), cex = cex.rates, srt = 90)
    text(x = px7+tdx, y = (py7+py8)/2.0, ratetext(x[7, 8], z[7, 8]), cex = cex.rates, srt = 90)
    text(x = px4-tdx, y = (py4+py6)/2.0, ratetext(x[6, 4], z[6, 4]), cex = cex.rates, srt = 90)
    text(x = px4+tdx, y = (py4+py6)/2.0, ratetext(x[4, 6], z[4, 6]), cex = cex.rates, srt = 90)

    text(x = (px1+px3)/2.0-tdx*0.7-amx, y = (py1+py3)/2.0+tdy*0.7+amy,
         ratetext(x[1, 3], z[1, 3]), cex = cex.rates, srt = 45)
    text(x = (px1+px3)/2.0+tdx*0.7-amx, y = (py1+py3)/2.0-tdy*0.7+amy,
         ratetext(x[3, 1], z[3, 1]), cex = cex.rates, srt = 45)
    text(x = (px4+px7)/2.0-tdx*0.7+amx, y = (py4+py7)/2.0-tdy*0.7+amy,
         ratetext(x[4, 7], z[4, 7]), cex = cex.rates, srt = -45)
    text(x = (px4+px7)/2.0+tdx*0.7+amx, y = (py4+py7)/2.0+tdy*0.7+amy,
         ratetext(x[7, 4], z[7, 4]), cex = cex.rates, srt = -45)
    text(x = (px2+px5)/2.0+tdx*0.7-amx, y = (py2+py5)/2.0+tdy*0.7-amy,
         ratetext(x[2, 5], z[2, 5]), cex = cex.rates, srt = -45)
    text(x = (px2+px5)/2.0-tdx*0.7-amx, y = (py2+py5)/2.0-tdy*0.7-amy,
         ratetext(x[5, 2], z[5, 2]), cex = cex.rates, srt = -45)
    text(x = (px6+px8)/2.0+tdx*0.7+amx, y = (py6+py8)/2.0-tdy*0.7-amy,
         ratetext(x[6, 8], z[6, 8]), cex = cex.rates, srt = 45)
    text(x = (px6+px8)/2.0-tdx*0.7+amx, y = (py6+py8)/2.0+tdy*0.7-amy,
         ratetext(x[8, 6], z[8, 6]), cex = cex.rates, srt = 45)

}

plotRatesZInd <- function(x, z, ...) {
    x[is.na(x)] <- 0
    # print(x)
    z[is.na(z)] <- 0
    # print(z)

    if (hasArg(signif)) 
        signif <- list(...)$signif
    else signif <- 3
    if (hasArg(main)) 
        main <- list(...)$main
    else main <- NULL
    if (hasArg(cex.main)) 
        cex.main <- list(...)$cex.main
    else cex.main <- 1.2
    if (hasArg(sub)) 
        sub <- list(...)$sub
    else sub <- NULL
    if (hasArg(cex.sub)) 
        cex.sub <- list(...)$cex.sub
    else cex.sub <- 1
    if (hasArg(cex.traits)) 
        cex.traits <- list(...)$cex.traits
    else cex.traits <- 0.9
    if (hasArg(cex.rates)) 
        cex.rates <- list(...)$cex.rates
    else cex.rates <- 0.8
    if (hasArg(lwd.by.rate))
        lwd.by.rate <- list(...)$lwd.by.rate
    else lwd.by.rate <- FALSE
    if (lwd.by.rate) {
        rates <- c(x[x > 0])
        LWD <- round(x/min(rates))
        if (hasArg(max.lwd)) 
            max.lwd <- list(...)$max.lwd
        else max.lwd <- 10
        LWD[LWD > max.lwd] <- max.lwd
    }
    else LWD <- matrix(2, nrow(x), ncol(x))
    
    options(repr.plot.width = 8, repr.plot.height = 4)
    plot.new()
    par(mar = c(1.1, 2.1, 3.1, 2.1))
    plot.window(xlim = c(0, 2), ylim = c(0, 1), asp = 1)
    mtext(
        if (!is.null(main))
            main
        else
            "Rates"
        , side = 3, adj = 0, line = 1.2,
        cex = cex.main)
    text(x = 0.0, y = 1.0, adj = 0,
        if (!is.null(sub))
            sub
        else
            "(X, Y, Z):"
        , cex = cex.sub)
    
    px <- c(0.4, 1.0, 1.6)
    py <- c(0.8, 0.2)
    adx <- 0.15
    ady <- 0.1
    adx2 <- 0.1
    ady2 <- 0.1
    amx <- 0.05
    amy <- 0.05
    aix <- 0.025
    aiy <- 0.025
    tdx <- 0.08
    tdy <- 0.08
    
    px1 <- px[1]
    px2 <- px[1]
    px3 <- px[2]
    px4 <- px[2]
    px5 <- px[3]
    px6 <- px[3]
    
    py1 <- py[1]
    py2 <- py[2]
    py3 <- py[1]
    py4 <- py[2]
    py5 <- py[1]
    py6 <- py[2]
        

    # Vertical
    
    # (0,*,*)(1347) <-> (1,*,*)(2568)
    arrows(x0 = px2-aix, y0 = py2+ady, x1 = px1-aix, y1 = py1-ady,
           lwd = max(LWD[2, 1], 1),
           lty = if (LWD[2, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px1+aix, y0 = py1-ady, x1 = px2+aix, y1 = py2+ady,
           lwd = max(LWD[1, 2], 1),
           lty = if (LWD[1, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # (*,0,*)(1246) <-> (*,1,*)(3578)
    arrows(x0 = px4-aix, y0 = py4+ady, x1 = px3-aix, y1 = py3-ady,
           lwd = max(LWD[3, 1], 1),
           lty = if (LWD[3, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px3+aix, y0 = py3-ady, x1 = px4+aix, y1 = py4+ady,
           lwd = max(LWD[1, 3], 1),
           lty = if (LWD[1, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # (*,*,0)(1235) <-> (*,*,1)(4678)
    arrows(x0 = px6-aix, y0 = py6+ady, x1 = px5-aix, y1 = py5-ady,
           lwd = max(LWD[4, 1], 1),
           lty = if (LWD[4, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px5+aix, y0 = py5-ady, x1 = px6+aix, y1 = py6+ady,
           lwd = max(LWD[1, 4], 1),
           lty = if (LWD[1, 4] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    
    ratetext <- function(x, z) {
        return(paste(round(x, signif), " (z = ", round(z, signif), ")", sep=""))
    }

    text(x = px1, y = py1, "(0,*,*)", cex = cex.traits)
    text(x = px2, y = py2, "(1,*,*)", cex = cex.traits)
    text(x = px3, y = py3, "(*,0,*)", cex = cex.traits)
    text(x = px4, y = py4, "(*,1,*)", cex = cex.traits)
    text(x = px5, y = py5, "(*,*,0)", cex = cex.traits)
    text(x = px6, y = py6, "(*,*,1)", cex = cex.traits)

    text(x = px1-tdx, y = (py1+py2)/2.0, ratetext(x[2, 1], z[2, 1]), cex = cex.rates, srt = 90)
    text(x = px1+tdx, y = (py1+py2)/2.0, ratetext(x[1, 2], z[1, 2]), cex = cex.rates, srt = 90)
    text(x = px3-tdx, y = (py3+py4)/2.0, ratetext(x[3, 1], z[3, 1]), cex = cex.rates, srt = 90)
    text(x = px3+tdx, y = (py3+py4)/2.0, ratetext(x[1, 3], z[1, 3]), cex = cex.rates, srt = 90)
    text(x = px5-tdx, y = (py5+py6)/2.0, ratetext(x[4, 1], z[4, 1]), cex = cex.rates, srt = 90)
    text(x = px5+tdx, y = (py5+py6)/2.0, ratetext(x[1, 4], z[1, 4]), cex = cex.rates, srt = 90)

}

plotRatesZIndZ <- function(x, z, ...) {
    x[is.na(x)] <- 0
    # print(x)
    z[is.na(z)] <- 0
    # print(z)

    if (hasArg(signif)) 
        signif <- list(...)$signif
    else signif <- 3
    if (hasArg(main)) 
        main <- list(...)$main
    else main <- NULL
    if (hasArg(cex.main)) 
        cex.main <- list(...)$cex.main
    else cex.main <- 1.2
    if (hasArg(sub)) 
        sub <- list(...)$sub
    else sub <- NULL
    if (hasArg(cex.sub)) 
        cex.sub <- list(...)$cex.sub
    else cex.sub <- 1
    if (hasArg(cex.traits)) 
        cex.traits <- list(...)$cex.traits
    else cex.traits <- 0.9
    if (hasArg(cex.rates)) 
        cex.rates <- list(...)$cex.rates
    else cex.rates <- 0.8
    if (hasArg(lwd.by.rate))
        lwd.by.rate <- list(...)$lwd.by.rate
    else lwd.by.rate <- FALSE
    if (lwd.by.rate) {
        rates <- c(x[x > 0])
        LWD <- round(x/min(rates))
        if (hasArg(max.lwd)) 
            max.lwd <- list(...)$max.lwd
        else max.lwd <- 10
        LWD[LWD > max.lwd] <- max.lwd
    }
    else LWD <- matrix(2, nrow(x), ncol(x))
    
    options(repr.plot.width = 8, repr.plot.height = 4)
    plot.new()
    par(mar = c(1.1, 2.1, 3.1, 2.1))
    plot.window(xlim = c(0, 2), ylim = c(0, 1), asp = 1)
    mtext(
        if (!is.null(main))
            main
        else
            "Rates"
        , side = 3, adj = 0, line = 1.2,
        cex = cex.main)
    text(x = 0.0, y = 1.0, adj = 0,
        if (!is.null(sub))
            sub
        else
            "(X, Y, Z):"
        , cex = cex.sub)
    
    px <- c(0.2, 1.2, 1.8)
    py <- c(0.8, 0.2)
    adx <- 0.15
    ady <- 0.1
    adx2 <- 0.1
    ady2 <- 0.1
    amx <- 0.05
    amy <- 0.05
    aix <- 0.025
    aiy <- 0.025
    tdx <- 0.08
    tdy <- 0.08
    
    px1 <- px[1]  # (0,0,*)
    px2 <- px[1]  # (1,0,*)
    px3 <- px[2]  # (0,1,*)
    px4 <- px[2]  # (1,1,*)
    px5 <- px[3]  # (*,*,0)
    px6 <- px[3]  # (*,*,1)
    
    py1 <- py[1]
    py2 <- py[2]
    py3 <- py[1]
    py4 <- py[2]
    py5 <- py[1]
    py6 <- py[2]
        

    # Vertical
    
    # (0,0,*)(14) <-> (1,0,*)(26)
    arrows(x0 = px2-aix, y0 = py2+ady, x1 = px1-aix, y1 = py1-ady,
           lwd = max(LWD[2, 1], 1),
           lty = if (LWD[2, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px1+aix, y0 = py1-ady, x1 = px2+aix, y1 = py2+ady,
           lwd = max(LWD[1, 2], 1),
           lty = if (LWD[1, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # (0,1,*)(37) <-> (1,1,*)(58)
    arrows(x0 = px4-aix, y0 = py4+ady, x1 = px3-aix, y1 = py3-ady,
           lwd = max(LWD[5, 3], 1),
           lty = if (LWD[5, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px3+aix, y0 = py3-ady, x1 = px4+aix, y1 = py4+ady,
           lwd = max(LWD[3, 5], 1),
           lty = if (LWD[3, 5] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # (*,*,0)(1235) <-> (*,*,1)(4678)
    arrows(x0 = px6-aix, y0 = py6+ady, x1 = px5-aix, y1 = py5-ady,
           lwd = max(LWD[4, 1], 1),
           lty = if (LWD[4, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px5+aix, y0 = py5-ady, x1 = px6+aix, y1 = py6+ady,
           lwd = max(LWD[1, 4], 1),
           lty = if (LWD[1, 4] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # Horizontal
    
    # (0,0,*)(14) <-> (0,1,*)(37)
    arrows(x0 = px1+adx, y0 = py1+aiy, x1 = px3-adx, y1 = py3+aiy,
           lwd = max(LWD[1, 3], 1),
           lty = if (LWD[1, 3] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px3-adx, y0 = py3-aiy, x1 = px1+adx, y1 = py1-aiy,
           lwd = max(LWD[3, 1], 1),
           lty = if (LWD[3, 1] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    # (1,0,*)(26) <-> (1,1,*)(58)
    arrows(x0 = px2+adx, y0 = py2+aiy, x1 = px4-adx, y1 = py4+aiy,
           lwd = max(LWD[2, 5], 1),
           lty = if (LWD[2, 5] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)
    arrows(x0 = px4-adx, y0 = py4-aiy, x1 = px2+adx, y1 = py2-aiy,
           lwd = max(LWD[5, 2], 1),
           lty = if (LWD[5, 2] == 0)
               "dashed"
           else "solid",
           length = 0.15, lend = 3, angle = 20)

    
    ratetext <- function(x, z) {
        return(paste(round(x, signif), " (z = ", round(z, signif), ")", sep=""))
    }

    text(x = px1, y = py1, "(0,0,*)", cex = cex.traits)
    text(x = px2, y = py2, "(1,0,*)", cex = cex.traits)
    text(x = px3, y = py3, "(0,1,*)", cex = cex.traits)
    text(x = px4, y = py4, "(1,1,*)", cex = cex.traits)
    text(x = px5, y = py5, "(*,*,0)", cex = cex.traits)
    text(x = px6, y = py6, "(*,*,1)", cex = cex.traits)

    text(x = px1-tdx, y = (py1+py2)/2.0, ratetext(x[2, 1], z[2, 1]), cex = cex.rates, srt = 90)
    text(x = px1+tdx, y = (py1+py2)/2.0, ratetext(x[1, 2], z[1, 2]), cex = cex.rates, srt = 90)
    text(x = px3-tdx, y = (py3+py4)/2.0, ratetext(x[5, 3], z[5, 3]), cex = cex.rates, srt = 90)
    text(x = px3+tdx, y = (py3+py4)/2.0, ratetext(x[3, 5], z[3, 5]), cex = cex.rates, srt = 90)
    text(x = px5-tdx, y = (py5+py6)/2.0, ratetext(x[4, 1], z[4, 1]), cex = cex.rates, srt = 90)
    text(x = px5+tdx, y = (py5+py6)/2.0, ratetext(x[1, 4], z[1, 4]), cex = cex.rates, srt = 90)
    
    text(x = (px1+px3)/2.0, y = py1+tdy, ratetext(x[1, 3], z[1, 3]), cex = cex.rates)
    text(x = (px1+px3)/2.0, y = py1-tdy, ratetext(x[3, 1], z[3, 1]), cex = cex.rates)
    text(x = (px2+px4)/2.0, y = py2+tdy, ratetext(x[2, 5], z[2, 5]), cex = cex.rates)
    text(x = (px2+px4)/2.0, y = py2-tdy, ratetext(x[5, 2], z[5, 2]), cex = cex.rates)

}
