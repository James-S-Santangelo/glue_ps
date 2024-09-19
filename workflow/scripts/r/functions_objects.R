my_theme <- theme_classic() +
            theme(axis.text = element_text(size = 15),
                  axis.title = element_text(size = 17),
                  legend.title = element_text(size = 15),
                  legend.text = element_text(size = 13)) +
            theme(panel.background = element_rect(fill = "transparent", colour = NA_character_), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.background = element_rect(fill = "transparent", colour = NA_character_),
                  legend.background = element_rect(fill = "transparent"),
                  legend.box.background = element_rect(fill = "transparent"),
                  legend.key = element_rect(fill = "transparent"))
