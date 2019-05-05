library(hexSticker)
img <- "logos/A_gran.png"
sticker(subplot=img, s_width=0.5, s_x=1, s_y=0.7,
        package="hSDM", p_size=30, p_x=1, p_y=1.4,
        h_fill=grey(0.7), h_color=grey(0.4),
        url="https://ecology.ghislainv.fr/hSDM",
        u_color="black",
        u_size=3,
        filename="man/figures/logo.png")