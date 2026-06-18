# Home page description

home_description <- tabPanel(title = "Home",
                             includeMarkdown("./www/Ratlas_home.md"),
                             tags$br(),
                             tags$div(
                               style = "text-align:center;",
                               tags$a(
                                 href   = "https://mapmyvisitors.com/web/1c5ld",
                                 title  = "Visit tracker",
                                 target = "_blank",
                                 rel    = "noopener noreferrer",
                                 tags$img(
                                   src = "https://mapmyvisitors.com/map.png?cl=080808&w=a&t=n&d=XUwzaDUmgBeEdg5Jxn4tF0dAaeAvpf8FszNMlO6XRBA&co=ffffff&ct=808080",
                                   alt = "Map of Ratlas visitors",
                                   style = "border:0;"
                                 )
                               )
                             ),
                             tags$br()
)
