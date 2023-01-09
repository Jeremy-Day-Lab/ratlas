# Home page description

home_description <- tabPanel(title = "Home",
                             includeMarkdown("./www/Ratlas_home.md"),
                             tags$br(),
                             HTML(
                               "<center>
      <script type='text/javascript' id='clustrmaps' src='//cdn.clustrmaps.com/map_v2.js?cl=ffffff&w=300&t=n&d=iIPr-XExkHg2hYOF1lT3nHdwwvsClutCG8OrW8Kfn9E&co=363636&cmo=d55e00&cmn=009e73'></script>
      </center>"
                             ),
      tags$br())
