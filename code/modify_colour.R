modify_colour <- function(hex, offset) {
    lookup <- data.frame(code = c(0:9, LETTERS[1:6]), 
                         int = 0:15)
    
    # hex_to_RGB
    red <- substr(hex, 2, 3)
    green <- substr(hex, 4, 5)
    blue <- substr(hex, 6, 7)
    
    codes <- sapply(2:7, function(x) substr(hex, x, x))
    
    RGB_val <- lookup$int[match(codes, lookup$code)] * rep(c(16,1), 3)
    RGB_val2 <- c(sum(RGB_val[1:2]), sum(RGB_val[3:4]), sum(RGB_val[5:6]))

    # apply shade/tint
    if(offset < 0) {
        RGB_val3 <- RGB_val2 * (1 - abs(offset))
    } 
    if(offset > 0) {
        RGB_val3 <- RGB_val2 + (255 - RGB_val2) * offset
    }
    if(offset == 0) {
        RGB_val3 <- RGB_val2
    }
    RGB_val3 <- round(RGB_val3, 0)
    
    # RGB to hex
    second_digit <- RGB_val3 %% 16
    first_digit <- (RGB_val3 - second_digit)/16
    
    hex_out <- c()
    hex_out[c(1, 3, 5)] <- first_digit
    hex_out[c(2, 4, 6)] <- second_digit
    
    paste0("#", paste0(lookup$code[match(hex_out, lookup$int)], collapse=""))
}
