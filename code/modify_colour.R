modify_colour <- function(hex, offset) {
    # hex <- "#FF0000"
    # hex <- c("#FF0000", "#00FF00")
    # offset <- c(-.5, .5)
    if(length(hex) == 1 & length(offset) >1) {
        hex <- rep(hex, length(offset))
    }
    if(length(hex) > 1 & length(offset) == 1) offset <- rep(offset, length(hex))

    if(length(hex) != length(offset)) stop("dimension mismatch")
    
    lookup <- data.frame(code = c(0:9, LETTERS[1:6]), 
                         int = 0:15)
    
    hex_mat <- do.call(rbind, strsplit(hex, ""))[,-1]
    
    
    # hex_to_RGB
    red <- hex_mat[,1:2]
    green <- hex_mat[,3:4]
    blue <- hex_mat[,5:6]
    
    codes <- sapply(2:7, function(x) substr(hex, x, x))

    # empty matrix
    RGB_val <- matrix(nrow=nrow(hex_mat), ncol=ncol(hex_mat))
    # populate and multiply
    RGB_val[] <- lookup$int[match(hex_mat, lookup$code)] 
    RGB_val <- RGB_val %*% diag(rep(c(16,1), 3))
    
    RGB_val2 <- data.frame(R = rowSums(RGB_val[,1:2]), 
                           G = rowSums(RGB_val[,3:4]), 
                           B = rowSums(RGB_val[,5:6]))
    
    # apply shade/tint
    index_shade <- offset < 0
    index_tint<- offset > 0
    
    RGB_val3 <- RGB_val2
    RGB_val3[index_shade,] <- RGB_val2[index_shade,] * (1 - abs(offset[index_shade]))
    RGB_val3[index_tint,] <- RGB_val2[index_tint,] + 
        (255 - RGB_val2[index_tint,]) * offset[index_tint]

    RGB_val3 <- round(RGB_val3, 0)
    
    # RGB to hex
    second_digit <- RGB_val3 %% 16
    first_digit <- (RGB_val3 - second_digit)/16
    
    hex_out <- cbind(first_digit[,1], second_digit[,1], 
                     first_digit[,2], second_digit[,2], 
                     first_digit[,3], second_digit[,3])
    
    lookup$code[match(hex_out, lookup$int)]
    hex_out[] <- lookup$code[match(hex_out, lookup$int)]
 
    paste0("#", do.call(paste0, data.frame(hex_out)))
}
