# this app is based on compare_gene_list.R

# remove user library path to avoid confusion
if(length(.libPaths())>1) .libPaths(.libPaths()[-1])
run_nextflow <- FALSE
library(shiny)
library(GeneOverlap)
library(gridExtra)
library(ggplot2)
library(gridGraphics)
library(eulerr)
library(reshape)
library(ggthemes)
library(ggpointdensity)
library(shiny)
library(RColorBrewer)
library(grid)
library(gridBase)
library(venn) # need to be 1.11
library(UpSetR)
library(DESeq2)
library(clusterProfiler) # need to be 3.8.1
library(msigdbr)
library(stringr)
library(ggplotify) #as.ggplot for upsetr
library(ggpolypath) #ggplot=true for venn
library(shinyFiles)
library(shinyjs)
library(ggrepel)
library(dplyr) #bind_rows, case_when
library(gtools) # mixedsort
#library(plotly)
#library(DT)
library(purrr)
library(shinyBS)
#library(ggigraph)
#library(tidyverse)

css <- "
.nav li a.disabled {
background-color: #d8d8d8 !important;
color: #8a8a8a !important;
cursor: not-allowed !important;
border-color: #aaa !important;
}"

#default hostfilepath to "/" if not exist
if(!exists("hostfilepath")) hostfilepath="/"
if(exists("hostfilepath") && substr(hostfilepath,nchar(hostfilepath),nchar(hostfilepath))!="/") 
 hostfilepath=paste0(hostfilepath,"/")

# default image directory used to run analysis out of container
#if(!exists("img.dir")) img.dir="~/Steve/virtual_server/rnaseq-singularity"

# ad-hoc variable to prevent error in initialization
log.lst=""
log.fail.lst=""


 # default checkboxes for venn.opts	
 venn.opts=data.frame(
  lbl= c("Numbers","Percentages","Labels","Legend")
 )
 venn.opts.lst=as.list(venn.opts$lbl)
 names(venn.opts.lst)=venn.opts$lbl
 
 # list of choices for msigdb species
 msigdb.species.lst <- as.list(msigdbr_species())$species_name
 
 ngroup=3
 tab0.ngroup=2
 
 colors=c("#F08080","#ADD8E6","#FFFACD")
 tab0.colors=c("#F08080","#ADD8E6","#FFFACD")
 tab3.colors=brewer.pal(12,"Paired")
 # color for pvalues (<2.2e-16,<=0.05,<0.5,>0.5)
 #col.pval=c("#fc8d59","#ffffbf","#91bfdb")
 col.pval = brewer.pal(11,"RdBu")[c(3,4,8,9)]
 # names for factorized pvalues
 pval.names=c("< 2.2e-16","< 0.05","< 0.5","> 0.5")
 # ad-hoc maximum plot. 
 max_plots=5
 # color pallette for pathways p-values
 pal <- colorRampPalette(c(tab3.colors[1],tab3.colors[2]))

# arbitrary variables to prevent errors
# groups.lst=c("a","b","c")
 
# List of choices for genome
setup.genome.lst <- as.list(c("hg19","hg38","other"))
# ad-hoc max.nsamples if not exist set
if(!exists("max.nsamples")) max.nsamples=50

heatmap_sigf_overlap <- function(data,title=""){
 hm <- ggplot(data,aes(X1,X2,fill=pval.cat)) +
  # use thin border size 0.1 to separate tiles
  geom_tile(color="white",size=0.1)+
  # setting high and low color and color key name
  scale_fill_manual(values=col.pval,na.value = "transparent",drop=FALSE)+
  # draw nice square
  # coord_equal() +
  # remove x and y labels
  labs(x=NULL, y=NULL,fill="p-values") +
  
  scale_x_discrete(labels = function(x) str_wrap(gsub("_"," ",x), width = 10)) +
  # removes a lot of chart junks including space around heatmap
  # from ggthemes package
  theme_tufte(base_family="Helvetica") +
  # remove tick marks on axes
  theme(axis.ticks=element_blank()) +
  # make text bigger
  theme(axis.text = element_text(size=14)) +
  # rotate x axis text
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # make legend on top left
  theme(legend.position="top") +
  #theme(legend.justification = c(-1,0)) +
  # add pvalue as text
  geom_text(aes(label=round(jaccard,digits=2))) +
  # color key size
  theme(legend.key.width=unit(0.7,"cm")) +
  ggtitle(title)
 # plot venn and overlap heatmap together
 grid.arrange(hm,
              top=textGrob(""),
              bottom=textGrob(paste0("Color represents p-values of overlap.\n",
                                     "Number inside a cell is Jaccard similarity index.\n",
                                     "0 means no similarity, 1 means identical."),just=c("left"),
                              hjust=0,x=0.3,gp=gpar(fontsize=12)))	 
 # grob <- grobTree(textGrob(paste("Color represents p-values of overlap.\n",
 # "Number inside a cell is Jaccard similarity index.\n",
 # "0 means no similarity, 1 means identical."), x=0,  y=0, hjust=0,
 # gp=gpar(fontsize=13, fontface="italic")))
 # hm + annotation_custom(grob)
}

plot_euler <- function(s4,colors,cex,venn.opts,title,labels){
 quantities=list(cex=cex);
 legend=FALSE
 if('Numbers' %in% venn.opts) quantities$type=c(quantities$type,"counts")
 if('Percentages' %in% venn.opts) quantities$type=c(quantities$type,"percent")
 if('Legend' %in% venn.opts) legend=list(side="bottom")
 
 if(is.null(venn.opts) | 
    sum(c('Numbers','Percentages') %in% venn.opts)==0) quantities=FALSE
 
 set.seed(1)
 plot(euler(s4,shape="ellipse"),
      quantities=quantities, labels=labels,
      fills=colors,edges=FALSE,main=title,
      adjust_labels = TRUE,legend=legend)
}

# Function to calculate directory size and return in human-readable format
get_dir_size <- function(path) {
 if (dir.exists(path)) {
  size <- sum(file.info(list.files(path, full.names = TRUE, recursive = TRUE))$size)
  # Convert size to human-readable format
  if (size < 1024) {
   size_str <- paste(size, "bytes")
  } else if (size < 1024^2) {
   size_str <- paste(round(size / 1024, 2), "KB")
  } else if (size < 1024^3) {
   size_str <- paste(round(size / 1024^2, 2), "MB")
  } else {
   size_str <- paste(round(size / 1024^3, 2), "GB")
  }
  return(size_str)
 } else {
  return("Directory does not exist.")
 }
}

ui <- fluidPage(
 title = "Run full RNA-Seq analysis",
 shinyjs::useShinyjs(),
 shinyjs::inlineCSS(css),
 # horizontal vertical bar for r1 and r2 filepaths input text
 # tags to have horizontal scrolling bar
 # white-space: pre will make newline character works. without this, newline chr in text will be ignored
 tags$style(HTML("#setup_grp1_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tags$style(HTML("#setup_grp1_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tags$style(HTML("#setup_grp2_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tags$style(HTML("#setup_grp2_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}")),
 tabsetPanel(id="tabset",
  tabPanel("Run Analysis",id="run_analysis", fluid = TRUE,
           sidebarPanel(width=3,
                        shinyDirButton("setup_proj_dir", 
                                         label = "Select project folder",
                                         title="Please select a folder to run pipeline",
                                       multiple=FALSE),
           bsButton("proj-folder-info", label = "", 
                    icon = icon("info", lib = "font-awesome"), style = "default", 
                    size = "extra-small"), 
  bsPopover(
   id = "proj-folder-info",
   title = "More information",
   content = HTML(paste0(
    "Select a folder which will contain all the files and outputs generated by the pipeline. ",
    "Avoid space in the path."
   )),
   placement = "right",
   trigger = "hover",
   options = list(container = "body")
  ),
  verbatimTextOutput('setup.proj.dir.filepaths'),
  br(),
  #verbatimTextOutput("fileExists"),
  conditionalPanel(condition="output.fileExists && !input.setup_proj_dir == '' ",
                   actionButton("setup.load.samples",
                                label=list("Click to load existing samples.txt")),
                   bsButton("load-sample-info", label = "", 
                            icon = icon("info", lib = "font-awesome"), style = "default", 
                            size = "extra-small"), 
                   bsPopover(
                    id = "load-sample-info",
                    title = "More information",
                    content = HTML(paste0("The existing samples.txt must be created ",
                                     "using the same filepath.")),
                    placement = "right",
                    trigger = "hover",
                    options = list(container = "body")
                   ),
                   br(),
                   br(),
  ),
  
  selectInput(inputId = "setup.genome",label="Select genome",
              choices=setup.genome.lst,selected=1),
  
  conditionalPanel(
   condition= "input['setup.genome'] == 'other'",
   textInput(inputId="setup.genome.name", 
             label = list("Genome name",
                          bsButton("setup-genome-name-info", label = "", 
                                   icon = icon("info", lib = "font-awesome"), 
                                   style = "default", size = "extra-small")),value=""),
   bsPopover(
    id = "setup-genome-name-info",
    title = "More information",
    content = HTML(paste0(
     "Name of genome assembly, e.g., mm10, danRer11, etc.<br> ",
     "STAR_index will be searched and/or created in the same folder as the FASTA file below, ",
     "within a subfolder named according to the Genome name specified by user."
    )),
    placement = "right",
    trigger = "hover",
    options = list(container = "body")
   ),
   
   shinyFilesButton("setup_genome_fa", 
                    label = "Select genome fasta file",
                                           title="Please select genome fasta file",multiple=FALSE),
                         bsButton("setup-genome-fa-info", label = "", 
                                  icon = icon("info", lib = "font-awesome"), 
                                  style = "default", size = "extra-small"), 
                         bsPopover(
                          id = "setup-genome-fa-info",
                          title = "More information",
                          content = HTML(paste0(
                           "Path to genome fasta file. Gzipped file is not supported."
                           )),
                          placement = "right",
                          trigger = "hover",
                          options = list(container = "body")
                         ),
                         verbatimTextOutput('setup.genome.fa.filepaths'),
                         br(),
                         shinyFilesButton("setup_genome_gtf", 
                                          label = "Select genome GTF file",
                                          title="Please select genome GTF file",multiple=FALSE ),
                         bsButton("setup-genome-gtf-info", label = "", icon = icon("info", lib = "font-awesome"), style = "default", size = "extra-small"), 
                         bsPopover(
                          id = "setup-genome-gtf-info",
                          title = "More information",
                          content = HTML(paste0(
                           "Path to genome annotation GTF file. Gzipped file is not supported.",
                           "Make sure annotation contains gene symbols. RefSeq annotation is preferred."
                          )),
                          placement = "right",
                          trigger = "hover",
                          options = list(container = "body")
                         ),
                         verbatimTextOutput('setup.genome.gtf.filepaths'),
                         br()
                        ),#conditionPanel
                        textInput(inputId = "setup.time",
                                  label=list("Time limit",
                                      bsButton("time-limit-info", label = "", icon = icon("info", lib = "font-awesome"), 
                                               style = "default", size = "extra-small")), 
                                                         value="7-00:00:00"),
  
  bsPopover(
   id = "time-limit-info",
   title = "More information",
   content = HTML(paste0(
    "Specify time limit DD-HH:MM:SS (D-days, H-hour, M-minutes, S-seconds)"
   )),
   placement = "right",
   trigger = "hover",
   options = list(container = "body")
  ),
                        checkboxInput(inputId = "setup.batch.adjust",
                                      label=list("Replicates batch adjustment",
                                              bsButton("setup-batch-info", label = "", 
                                              icon = icon("info",lib = "font-awesome"),
                                              style = "default", size = "extra-small")), 
                                      value=TRUE),
  bsPopover(
   id = "setup-batch-info",
   title = "More information",
   content = HTML(paste0(
    "Batch correction to account for unwanted variation across different replicates.<br>",
    "Uncheck if number of replicates is not the same across groups."
   )),
   placement = "right",
   trigger = "hover",
   options = list(container = "body")
  ),
                        numericInput(inputId = "setup.ncpus.trim",
                                     label=list("# of CPUs for trimming",
                                                bsButton("setup-ncpus-trim-info", label = "", 
                                                icon = icon("info", lib = "font-awesome"), 
                                                style = "default", size = "extra-small")),
                                     value=4,min=1),
  bsPopover(
   id = "setup-ncpus-trim-info",
   title = "More information",
   content = HTML(paste0(
    "Number of CPUs to use for trimming fastq files.<br>",
    "Set to 1 if pigz error occurs."
   )),
   placement = "right",
   trigger = "hover",
   options = list(container = "body")
  ),
                        numericInput(inputId = "setup.ncpus.star",
                                     label=list("# of CPUs for STAR alignment",
                                                bsButton("setup-ncpus-star-info", label = "", 
                                                         icon = icon("info", lib = "font-awesome"), 
                                                         style = "default", size = "extra-small")), 
                                     value=20,min=1),
  bsPopover(
   id = "setup-ncpus-star-info",
   title = "More information",
   content = HTML(paste0(
    "Number of CPUs to use for STAR alignment."
   )),
   placement = "right",
   trigger = "hover",
   options = list(container = "body")
  ),
                        checkboxInput(inputId = "setup.debug",
                                      label=list("Run debug",
                                                 bsButton("setup-debug-info", label = "", 
                                                                      icon = icon("info", lib = "font-awesome"), 
                                                                      style = "default", size = "extra-small")),
                                      value=FALSE),
  bsPopover(
   id = "setup-debug-info",
   title = "More information",
   content = HTML(paste0(
    "Run with more messages printed to help with debugging errors."
   )),
   placement = "right",
   trigger = "hover",
   options = list(container = "body")
  ),
  actionButton(inputId="setup.save.samples","Save samples.txt"),
  actionButton(inputId="setup.view.samples","View samples.txt"),
  br(),
  br(),
                        conditionalPanel(condition = "input.setup.proj.dir == ''",
                        textOutput('err_msg')),
                        #conditionalPanel(condition = "input.setup.genome == 'other' && (input.setup.genome.name == '' || input.setup_genome_fa == '' || input.setup_genome_gtf == '')",
                                         textOutput('err_msg_genome'),
                        actionButton(inputId="setup.run.analysis",
                                     label=list("Run full analysis"),
                        ),
  bsButton("setup-run-analysis-info", label = "", icon = icon("info", lib = "font-awesome"), 
           style = "default", size = "extra-small"), 
  
  bsPopover(
   id = "setup-run-analysis-info",
   title = "More information",
   content = HTML(paste0(
    "Run full RNA-seq analysis: QC, trimming, alignment, tracks creation, reads counting, and differential gene analysis.<br>",
    "Once clicked, you must restart from the command line to perform a new analysis."
   )),
   placement = "right",
   trigger = "hover",
   options = list(container = "body")
  ),
  shinyjs::hidden(p(id="text1","Running full analysis. See Log tab for more info."))
           ),#sidebarPanel
           mainPanel(
            # adding horizontal scrolling
            style = "overflow-x: auto;", # Apply CSS styles
            width = 9,
                        fluidRow(
                         h3("Creating a samples.txt"),
                         column(4, 
                          textInput(inputId="setup.email", 
                                  label = list("Email address",
                                               bsButton("setup-email-info",label="",icon = icon("info", lib = "font-awesome"),
                                                        style = "default", size = "extra-small")), value=""),
                                  )
                                         ,
                         bsPopover(
                          id = "setup-email-info",
                          title = "More information",
                          content = HTML(paste0(
                           "Enter your email address to get emails for completed or failed jobs.")),
                          placement = "right",
                          trigger = "hover",
                          options = list(container = "body")
                         )
                         ),
                     fluidRow(
                      column(2,tags$label("Group name",
                                          bsButton("setup-group-info", label = "", icon = icon("info", lib = "font-awesome"), 
                                                   style = "default", size = "extra-small"))
                      ,style="padding-left:30px;padding-top:0px"),
                      bsPopover(
                       id = "setup-group-info",
                       title = "More information",
                       content = HTML(paste0(
                        "Group name of the samples. Samples with replicates must have the same Group name. ",
                        "This column cannot be empty.")),
                       placement = "right",
                       trigger = "hover",
                       options = list(container = "body")
                      ),
                      column(2,tags$label("Control name",
                       bsButton("setup-control-info", label = "", icon = icon("info", lib = "font-awesome"), style = "default", size = "extra-small"))
                      ),
                      bsPopover(
                       id = "setup-control-info",
                       title = "More information",
                       content = HTML(paste0(
                        "Name of the control/background/reference sample to be compared to group name sample. ",
                        "During differential analysis, ",
                        "the corresponding sample will be compared to its Control name. ",
                        "Fold-changes will be calculated with Control name as reference. ",
                        "Control name cannot be empty. For control/background sample, please put NA"
                       )),
                       placement = "right",
                       trigger = "hover",
                       options = list(container = "body")
                      ),
                      column(2,tags$label("Replicate name",
                           bsButton("setup-rep-info", label = "", icon = icon("info", lib = "font-awesome"), 
                                    style = "default", size = "extra-small"))
                      ),
                      bsPopover(
                       id = "setup-rep-info",
                       title = "More information",
                       content = HTML(paste0(
                        "Name of the replicates.",
                        "Avoid naming rep1 as R1 as that is reserved as read pair 1. ",
                        "Group name, followed by rep name will be your output filenames."
                       )),
                       placement = "right",
                       trigger = "hover",
                       options = list(container = "body")
                      ),
                      column(3,tags$label("Path to R1 fastq files",
                                          bsButton("setup-r1-path-info", 
                                                   label = "", icon = icon("info", lib = "font-awesome"), 
                                                   style = "default", size = "extra-small"))
                      ),
                      bsPopover(
                       id = "setup-r1-path-info",
                       title = "More information",
                       content = HTML(paste0(
                        "Select the fastq file for read pair 1. ",
                        "To select multiple files, hold ctrl button in windows. ",
                        "Click the second time to select additional files in different path. ",
                        "You can also change the path manually below or copy paste paths. ",
                        "Do not use path with white spaces."
                       )),
                       placement = "right",
                       trigger = "hover",
                       options = list(container = "body")
                      ),
                      column(3,tags$label("Path to R2 fastq files",
                                          bsButton("setup-r2-path-info", 
                                                   label = "", icon = icon("info", lib = "font-awesome"), 
                                                   style = "default", size = "extra-small"))
                     ),
                     bsPopover(
                      id = "setup-r2-path-info",
                      title = "More information",
                      content = HTML(paste0(
                       "Select the fastq file for read pair 2. ",
                       "To select multiple files, hold ctrl button in windows. ",
                       "Click the second time to select additional files in different path. ",
                       "You can also change the path manually below or copy paste paths. ",
                       "Do not use path with white spaces."
                      )),
                      placement = "right",
                      trigger = "hover",
                      options = list(container = "body")
                     )),
                     hr(style="border-color:black;margin-top:0px;margin-bottom:0px"),
                     fluidRow(style = "background-color:#f9f9f9",
                              column(2,
                                     column(2, 
                                            tags$label("1.", style = "padding-top: 30px;"),
                                            style = "padding-right: 0px;padding-left:0px"),
                                     column(10,
                                 
                                           textAreaInput(inputId="setup.grp1.name", 
                                                         label = "",value=""),
                                           style ="padding-right:0px;padding-left:0px")
                              ),
                              column(2,textAreaInput(inputId="setup.grp1.ctrl.name", 
                                                 label = "",value="" )),
                              column(2,textInput(inputId="setup.grp1.rep.name", 
                                                 label = "",value="" )),
                              br(),
                              column(3,shinyFilesButton("setup_grp1_r1_files", 
                                                        label = "Select R1 fastq files",
                                                        title="Please select Group 1 R1 fastq files",multiple=TRUE ),
                                     textAreaInput('setup_grp1_r1_filepaths',label="",value="",width='150px',height='100px')),
                              column(3,shinyFilesButton("setup_grp1_r2_files", 
                                                        label = "Select R2 fastq files",
                                                        title="Please select Group 1 R2 fastq files",multiple=TRUE ),
                                     textAreaInput('setup_grp1_r2_filepaths',label="",width='150px',height='100px')),
                     ), #fluidRow
                     fluidRow(
                      column(2,column(2, 
                                      tags$label("2.", style = "padding-top: 30px;"),
                                      style = "padding-right: 0px;padding-left:0px"),
                             column(10,textAreaInput(inputId="setup.grp2.name", 
                                                 label = "",value="" ),
                                    style ="padding-right:0px;padding-left:0px")
                      ),
                      column(2,textAreaInput(inputId="setup.grp2.ctrl.name", 
                                         label = "",value="" )),
                      column(2,textInput(inputId="setup.grp2.rep.name", 
                                         label = "",value="" )),
                      br(),
                      column(3,shinyFilesButton("setup_grp2_r1_files", 
                                                label = "Select R1 fastq files",
                                                title="Please select Group 2 R1 fastq files",multiple=TRUE ),
                             textAreaInput('setup_grp2_r1_filepaths',label="",value="",width='150px',height='100px')),
                                 column(3,shinyFilesButton("setup_grp2_r2_files", 
                                                label = "Select R2 fastq files",
                                                title="Please select Group 2 R2 fastq files",multiple=TRUE ),
                                        textAreaInput('setup_grp2_r2_filepaths',label="",value="",width='150px',height='100px')),
                     ), #fluidRow
                     # tag$div id doesn't work if there is a dot!
                     tags$div(id = "setup-placeholder"),
                     hr(style="border-color:black;margin-top:0px;margin-bottom:0px"),
                     br(),
                     column(6,actionButton("setup.remove.set", "Remove row", width = "100%"),
                            style="padding-bottom:20px"),
                     column(6,actionButton("setup.insert.set", "Add row", width = "100%"),
                            style="padding-bottom:20px"),
           ),
           
  ),# tabPanel run analysis
  tabPanel("Log", id="Log", fluid = TRUE,
           sidebarPanel(width=3,
                        # Input: Choose a log file to view
                        actionButton("logtab.refresh.log.path","Refresh list"),
                        br(),
                        br(),
                        selectInput("logtab.log.path", "Choose a log file to view:",
                                    choices = log.lst,selected="run_rnaseq_full.out"),
                        # commenting this out since I can't find FAILED word in other slurm settings
                        #br(),
                        ## Input: Choose a failed log file to view
                        #selectInput("logtab.fail.log.path", "Choose a failed log file to view:",
                        #            choices = log.fail.lst,selected="")
                        ),
           mainPanel(width=9,
                     fluidPage(
                      tags$div(
                      h5("Log file content:"),
                      tags$head(
                       tags$style(HTML(
                       "#logtab_log_content{ /* id of verbatimTextOutput CAN'T have dot */
                        height: 90vh;  /* Set height to 90% of the viewport height */
                        overflow-y: auto;  /* Enable vertical scrolling */
                        overflow-x: auto;  /* Enable horizontal scrolling */
                       }"
                       ))
                       ),
                      verbatimTextOutput("logtab_log_content"),
                     )#div
                     # commenting this out since I can't find FAILED word in other slurm settings
                      #h5("Failed log content:"),
                     #verbatimTextOutput("logtab.fail.log.content")
                     ))
                     
  ),# tabPanel Log
  tabPanel("QCs",id="qctab",fluid=TRUE,
            #uiOutput("multiqc")
           htmlOutput("multiqc")
  ), #tabPanel QCs
  tabPanel("Outputs",id="outputtab",fluid=TRUE,
           fluidPage(
            htmlOutput("diff_report")
           )#fluidPage
            ), #tabPanel Outputs
  tabPanel("Plots", id="tab3",fluid = TRUE,
           sidebarLayout(
            sidebarPanel(
             tags$style(type='text/css', 
                        ".selectize-input { 
					 font-size: 10pt; 
					 line-height: 10pt;} 
					 .selectize-dropdown {
					 font-size: 10pt; 
					 line-height: 10pt;}"),
             width=2,
             
             splitLayout(
              numericInput("pdf.width", label = HTML("PDF <br/> width (in)"), 
                           value = 11),
              numericInput("pdf.height", label = HTML("PDF <br/> height (in)"), 
                           value = 8)),
             actionButton(inputId='export', label=list("Save plots.pdf"),),
             bsButton("export-info", label = "", 
                      icon = icon("info", lib = "font-awesome"), 
                      style = "default", size = "extra-small"),
             bsPopover(
              id = "export-info",
              title = "More information",
              content = HTML(paste0(
               "Plots are saved in your &ltproject folder&gt/outputs/shiny_outputs."
              )),
              placement = "right",
              trigger = "hover",
              options = list(container = "body")
             ), 
             br(),
             br(),
             actionButton(inputId='exportGene', label=list("Export filtered gene list"),),
                                  bsButton("exportGene-info", label = "", 
                                  icon = icon("info", lib = "font-awesome"), 
                                  style = "default", size = "extra-small"),
             bsPopover(
              id = "exportGene-info",
              title = "More information",
              content = HTML(paste0(
               "Filtered gene list is saved in your &ltproject folder&gt/outputs/shiny_outputs."
              )),
              placement = "right",
              trigger = "hover",
              options = list(container = "body")
             ) 
            ),
            mainPanel(width=10,
             fluidRow(style = "background-color:#F5F5F5",
                      column(2,selectInput(
                       inputId="tab3.grp1.name",
                       label = ("Group 1"),
                       choices = NULL),
                       tags$style(    type = 'text/css',
                                      ".selectize-input { word-wrap : break-word;}
                         .selectize-dropdown {word-wrap : break-word;} "
                       )),
                      column(2,textAreaInput(inputId="tab3.grp1.plot.title", 
                                label = ("Group 1 label"), 
                                value = "")),
                     column(2,
                             #textInput("tab3.color.grp1", label = HTML("Color for <br/> group 1"), value = tab3.colors[1]),
                             textInput("tab3.color.grp1", label = ("Color for group 1"), value = tab3.colors[1])),
                      column(2,
                             numericInput("tab3.fdr.grp1", label = ("FDR cut-off for group 1"), 
                                          value = 0.05,step=0.01,min=0,max=1)),
                      column(2,
                             numericInput("tab3.fc.grp1", label = ("Fold-change for group 1"), 
                                          value = 1,min=1,step=0.5)),
                      column(2,
                             numericInput("tab3.meanDiff.grp1", label = ("Difference for group 1"), 
                                          value = 0,min=0,step=1))),
                      tags$div(id = "placeholder"),
             splitLayout(
              actionButton("insert_set", "Insert", width = "100%"),
              actionButton("remove_set", "Remove", width = "100%")
             ),
             br(),
             br(),
             tabsetPanel(
              tabPanel("Table(s)",
             # Table of diff genes
             fluidRow(style = "background-color:#F5F5F5",
              br(),
              p("Table of gene expression treatment_vs_reference. 
                 Up-regulated means the expression is higher in treatment compared to reference.
                 NS means not significant.",style="font-size:18px"),
              textInput("tab3.search.term", "Search by gene names separated by commas:"),
              numericInput("currentPage","Current Page",value=1,min=1),
              
              uiOutput("tables")
             )),
             

 
             tabPanel("Volcano Plot(s)",
             fluidRow(
              br(),
              # genes to highlight in volcano plot
              column(5,textInput("tab3.hilite.genes",
                                 label = "Enter official gene names separated by comma:",value="")),
               # This is the dynamic UI for volcano plots
               column(12,uiOutput("tab3.volcano.plots"))
              )),
              tabPanel("Overlaps",
                       br(),
                       column(2,
                              checkboxGroupInput("tab3.venn.opts", label = ("Display on Venns"),
                                                 choices = venn.opts.lst,
                                                 selected = venn.opts$lbl[1:3])),
                              
                              # slider bar font size
                       column(2,sliderInput("tab3.venn.cex", label = ("Venns' font size"), min = 0,
                                          max = 6, value = 1,step = 0.1)),
                              # slider bar pad
                       column(2,sliderInput("tab3.venn.pad", label = ("Venns' circle size"), min = 0,
                                          max = 8, value = 1,step = 0.5)),
                       br(),
               column(12,uiOutput("tab3.venn.plots"))
              ),
             tabPanel("Upset Plot(s)",
                      br(),
                             # upset plot number of intersection
                             numericInput("tab3.nintersects.upset", label = ("Upset plot max intersections"),
                                          value = 40,min=1,step=1),
                      
              plotOutput(outputId = "tab3.plot1")
             ,
             br(),
             br(),
              column(12,plotOutput(outputId = "tab3.plot2"))
              ),
             tabPanel("Pathway",
                      fluidRow(
                       br(),
             # species for msigdb
             selectInput(
              inputId="tab3.msigdb.species",
              label = ("Select your species"),
              selected="Homo sapiens",
              choices = msigdb.species.lst),
             # gene set enrichment pvalue numeric input
             numericInput("tab3.enrich.pval.co", label = ("FDR for enrichment"), 
                          value = 0.05,min=0,step=0.01),
             # max # of categories for pathway
             numericInput("tab3.enrich.ncat", 
                          label = list("# of categories",
                                       bsButton("tab3-enrich-ncat-info", label = "", 
                                                icon = icon("info", lib = "font-awesome"), 
                                                style = "default", size = "extra-small")), 
                          value = 10,min=1,step=1),
             bsPopover(
              id = "tab3-enrich-ncat-info",
              title = "More information",
              content = HTML(paste0(
               "Number of categories to show in Pathways. Note: the number of categories to show also depends on FDR."
              )),
              placement = "right",
              trigger = "hover",
              options = list(container = "body")
             ),
             textInput("tab3.chg.enrich.terms", label = ("Enter terms to change in enrichment plots (e.g., Gtpase=GTPase, Hiv=HIV)"),
                                 value=""),
             # Adjust pathways plot vertical size
             numericInput("tab3.enrich.plot.size", 
                          label = list("Plot vertical size",
                                       bsButton("tab3-enrich-plot-size-info", label = "", 
                                                icon = icon("info", lib = "font-awesome"), 
                                                style = "default", size = "extra-small")), 
                          value = 10,min=1,step=1),
             bsPopover(
              id = "tab3-enrich-plot-size-info",
              title = "More information",
              content = HTML(paste0(
               "Increase the number to increase vertical space in plots\n.",
               "Decrease to make it smaller."
              )),
              placement = "right",
              trigger = "hover",
              options = list(container = "body")
             ),
             actionButton(
              inputId = "gen.go",
              label = "Generate Enrichment plots"
             ),
             br(),
             br(),
             conditionalPanel(
              condition= "input['gen.go'] >= 1",
              fluidRow(
               column(12,plotOutput(outputId = "tab3.plot3.1",height="auto")),
         
              br(),
               column(12,plotOutput(outputId = "tab3.plot3.2",height="auto")),
              br(),
               column(12,plotOutput(outputId = "tab3.plot3.3",height="auto")),
              br(),
               column(12,plotOutput(outputId = "tab3.plot3.4",height="auto")),
              br(),
              column(12,plotOutput(outputId = "tab3.plot3.5",height="auto"))),
              ))))
            )
           )
  ), #tabPanel	plots
  tabPanel("Clean up", id="cleanuptab",fluid = TRUE,
           fluidPage(
            br(),
            shinyDirButton("del_folder", "Select project folder to clean up", "Please select a folder"),
            textOutput("selected_folder"),  # Display the selected folder path
            br(),
            div(
             style = "border-style: solid; border:1px solid rgba(0,0,0,0.15);
                padding:9.5px",
             p("Use \"Delete project folder\" to delete the entire folder"),
             actionButton("delete", "Delete project folder"),
             br(),
             textOutput("status")
            ),
            br(),
            div(
             style = "border-style: solid; border:1px solid rgba(0,0,0,0.15);
                padding:9.5px",
             p("Use \"Delete all files except samples.txt\" to delete all files inside
                      the project directory but keep samples.txt"),
             actionButton("delete_keep_samples", "Delete all files except samples.txt"),
             br(),
             textOutput("delete_keep_samples_status")
            ),#div
            br(),
            div(
             style = "border-style: solid; border:1px solid rgba(0,0,0,0.15);
                padding:9.5px",
             uiOutput("dir_size_merged_fastq"),
             p(HTML(paste(
               "Merged and trimmed fastq files were only used to re-run the same samples.<br>",
               "Typically, there is no need to save merged and trimmed fastq files unless you",
               "need to re-run the same exact analysis.<br>",
               "You should delete merged and trimmed fastq files once you're done 
                with the analysis.<br>"))
               ),
             p("Use \"Delete merged and/or trimmed fastq files\" to delete merged and trimmed fastq files inside
                      the project directory. All other files will not be deleted"),
             actionButton("delete_merged_fastq", "Delete merged and/or trimmed fastq files"),
             br(),
             textOutput("delete_merged_fastq_status")
            ),#div
            br(),
            div(
             style = "border-style: solid; border:1px solid rgba(0,0,0,0.15);
                padding:9.5px",
             uiOutput("dir_size_bam"),
             p(HTML(paste(
             "BAM files can be big. Keep your bam files if you are going to re-run the same",
             "samples. Otherwise, depending on your space limitation, you may or may not",
             "want to keep your bam files.<br>",
             "This pipeline generates bigwig files which are much smaller than bam files",
             "while still displaying read coverage.<br>",
             "Keep bam files if you have the storage capacity and want to be able to",
             "look at mapping of reads to the reference genome.<br>",
             "Otherwise, you should delete your bam files once you are done",
             "with your analysis.<br>",
             "Use \"Delete bam files\" to delete all bam files inside
                      the project directory"))),
             actionButton("delete_bam", "Delete all bam files"),
             br(),
             textOutput("delete_bam_status")
            ),#div
            br(),
           )#fluidpage
           )#tabpanel
  ,selected="Run Analysis")#tabsetpanel
)#fluidpage



# server
server <- function(input, output,session) {

 shinyjs::disable(selector = 'a[data-value="Two Groups"')
 shinyjs::disable(selector = 'a[data-value="Three Groups"')
 #shinyjs::disable(selector = 'a[data-value="Plots"')


 #write.table(isolate(session$clientData$url_port),file="/mnt/outputs/port.txt",
 #	quote=FALSE,row.names=FALSE,col.names=FALSE)
 
 #### run analysis tab ######
 # reactiveVal to keep track of number of added row
 # start with 2 as the minimum # of row
 setup.value <- reactiveVal(2)
 # initialize inserted div
 setup.inserted.div <- c()
 
 # inserting row
 observeEvent(input$setup.insert.set, {
  setup.btn <- setup.value() +1
  # update the reactiveVal
  setup.value(setup.btn)
  # id for new div
  setup.div.id <- paste0("setup_div_",setup.btn)
  # ids for entries in each row
  setup.id.grp.name <- paste0("setup.grp",setup.btn,".name")
  setup.id.ctrl.name <- paste0("setup.grp",setup.btn,".ctrl.name")
  setup.id.rep.name <- paste0("setup.grp",setup.btn,".rep.name")
  setup.id.r1.files <- paste0("setup_grp",setup.btn,"_r1_files")
  setup.id.r2.files <- paste0("setup_grp",setup.btn,"_r2_files")
  setup.id.r1.filepath <- paste0("setup_grp",setup.btn,"_r1_filepaths")
  setup.id.r2.filepath <- paste0("setup_grp",setup.btn,"_r2_filepaths")
  if(setup.btn %% 2 == 0) row.color="#FFFFFF"
  if(setup.btn %% 2 != 0) row.color="#f5f5f5"
  insertUI(
   selector = "#setup-placeholder",
   ui = tags$div(
    tags$style(HTML(paste0("#setup_grp",setup.btn,"_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}"))),
    tags$style(HTML(paste0("#setup_grp",setup.btn,"_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}"))),
    fluidRow(style = paste0("background-color:",row.color),
             column(2,
                    column(2,
                           tags$label(paste0(setup.btn,"."), style = "padding-top: 30px;"),
                           style = "padding-right: 0px;padding-left:0px"),
                    column(10,textAreaInput(inputId=setup.id.grp.name,
                                        label = "",value=""),
                                        style = "padding-right: 0px;padding-left:0px")
             ),
             column(2,textAreaInput(inputId=setup.id.ctrl.name,
                                label = "",value="" )),
             column(2,textInput(inputId=setup.id.rep.name,
                                label = "",value="" )),
             br(),
             column(3,shinyFilesButton(setup.id.r1.files,
                                       label = "Select R1 fastq files",
                                       title=paste("Please select Group",setup.btn,"R1 fastq files"),
                                       multiple=TRUE ),
                    textAreaInput(setup.id.r1.filepath,label="",value="",width='150px',height='100px')),
             column(3,shinyFilesButton(setup.id.r2.files,
                                       label = "Select R2 fastq files",
                                       title=paste("Please select Group",setup.btn,"R1 fastq files"),
                                       multiple=TRUE ),
             textAreaInput(setup.id.r2.filepath,label="",value="",width='150px',height='100px')),
    
    ),#fluidrow
    id = setup.div.id
   )
  )#insertUI
  setup.inserted.div <<- c(setup.inserted.div, setup.div.id)
 }) #observeEvent input$setup.insert.set
 
 observeEvent(input$setup.remove.set,{
  if(setup.value()>2){
   setup.btn <- setup.value()-1
   # update the reactiveVal
   setup.value(setup.btn)
   removeUI(
    selector = paste0("#",setup.inserted.div[length(setup.inserted.div)])
   )
   # remove the last values
   updateTextAreaInput(
    session,
    paste0("setup.grp", length(setup.inserted.div)+2,".name"),
    NULL,
    ""
   )
   updateTextAreaInput(
    session,
    paste0("setup.grp", length(setup.inserted.div)+2,".ctrl.name"),
    NULL,
    ""
   )
   updateTextInput(
    session,
    paste0("setup.grp", length(setup.inserted.div)+2,".rep.name"),
    NULL,
    ""
   )
   updateTextAreaInput(
    session,
    paste0("setup_grp", length(setup.inserted.div)+2,"_r1_files"),
    NULL,
    ""
   )
   updateTextAreaInput(
    session,
    paste0("setup_grp", length(setup.inserted.div)+2,"_r2_files"),
    NULL,
    ""
   )
   updateTextAreaInput(
    session,
    paste0("setup.grp", length(setup.inserted.div)+2,".r1.filepath"),
    NULL,
    ""
   )
   updateTextAreaInput(
    session,
    paste0("setup.grp", length(setup.inserted.div)+2,".r2.filepath"),
    NULL,
    ""
   )
   setup.inserted.div <<- setup.inserted.div[-length(setup.inserted.div)]
  }else{ # Do not remove row if there's only 1 row
   setup.btn <- 2
   setup.value(setup.btn)
  }
 })
 
  react.setup.grp.name <- reactive({
  grp.name <- sapply(grep("setup\\.grp\\d+\\.name", x = names(input), value = TRUE),
                     function(x) input[[x]])
  grp.name=as.vector(grp.name[mixedorder(names(grp.name))])
  grp.name[grp.name!=""]
 })
 
 react.setup.grp.ctrl.name <- reactive({
  grp.ctrl.name <- sapply(grep("setup\\.grp\\d+\\.ctrl\\.name", x = names(input), value = TRUE),
                          function(x) input[[x]])
  grp.ctrl.name=as.vector(grp.ctrl.name[mixedorder(names(grp.ctrl.name))])
  grp.ctrl.name[grp.ctrl.name!=""]
 })
 
 react.setup.grp.rep.name <- reactive({
  grp.rep.name <- sapply(grep("setup\\.grp\\d+\\.rep\\.name", x = names(input), value = TRUE),
                         function(x) input[[x]])
  # need to do order otherwise the new one will be the first in the vector
  grp.rep.name=as.vector(grp.rep.name[mixedorder(names(grp.rep.name))])
  grp.rep.name[grp.rep.name!=""]
 })
 
 react.setup.email <- reactive({
  input$setup.email
 })
 
 react.setup.genome.name <- reactive({
  input$setup.genome.name
 })
 
 react.setup.genome <- reactive({
  input$setup.genome
 })
 
 react.setup.genome.fa <- reactive({
  fa.file = input$setup_genome_fa
  fa.path=as.character(parseFilePaths(root = volumes,fa.file)$datapath)
  fa.path
 })
 
 react.setup.genome.gtf <- reactive({
  gtf.file = input$setup_genome_gtf
  gtf.path=as.character(parseFilePaths(root = volumes,gtf.file)$datapath)
  gtf.path
 })
 
 
 # browse from /filepath if it's specified otherwise /root
 if(file.exists("/filepath")){
  volumes=c(root="/filepath/")
 }else{ 
  volumes=c(root="/root/")
 }
 
 # shinyFileChoose doesn't work inside a reactive
 # I had to arbitrarily set a max number of samples a
 # samples.txt can have
 for (i in 1:max.nsamples) {
  setup.grp.r1.name=paste0('setup_grp',i,'_r1_files')
  shinyFileChoose(input, setup.grp.r1.name, root=volumes,
                  filetypes=c('', 'gz'),session=session)
  setup.grp.r2.name=paste0('setup_grp',i,'_r2_files')
  shinyFileChoose(input, setup.grp.r2.name, root=volumes,
                  filetypes=c('', 'gz'),session=session)
 }
 
 
 # Initialize r1/r2 file list as reactiveValues to save previous choices
 # need to use ad-hoc max.nsamples for shinyFile Choose
 x <- as.list(rep("",max.nsamples))
 names(x)=lapply(1:max.nsamples,function(i)paste0('grp',i))
 setup.grp.r1.file.lst <- do.call("reactiveValues",x)
 setup.grp.r2.file.lst <- do.call("reactiveValues",x)

  # getting all the r1 and r2 fastq files for all groups
  lapply(
  1:max.nsamples,
  function(i){
   observeEvent(input[[paste0('setup_grp',i,'_r1_files')]],{
    # getting new selected path and clean-up paths
    file = input[[paste0('setup_grp',i,'_r1_files')]]
    new.path=as.character(parseFilePaths(root = volumes,file)$datapath)
    new.path=gsub('^/filepath/','',new.path)
    new.path=gsub('^/root/','',new.path)
    # getting path in textAreaInput
    text.path=input[[paste0('setup_grp',i,'_r1_filepaths')]]
    text.path=unlist(strsplit(text.path,"\n"))
    # combine new.path and paths in textarea
    setup.grp.r1.file.lst[[paste0('grp',i)]] <- sort(unique(c(new.path,text.path)))
    path.display=gsub('^/filepath/','',setup.grp.r1.file.lst[[paste0('grp',i)]])
    print("path.display")
    print(path.display)
    path.display=gsub('^/root/','',path.display)
    print("path.display")
    print(path.display)
    # remove empty path
    if(sum(path.display=="")>0){
         rem.id=which(path.display=="")
         path.display=path.display[-rem.id]
    }
     
    path.display <- paste0(path.display, collapse="\n")
    updateTextAreaInput(session, paste0('setup_grp',i,'_r1_filepaths'),value=path.display)
    # output[[paste0('setup.grp',i,'.r1.filepaths')]]<- renderText({
    #  path.display=gsub('^/filepath/','',setup.grp.r1.file.lst[[paste0('grp',i)]])
    #  path.display=gsub('^/root/','',path.display)
    #  path.display <- paste0(path.display[-1], "\n")
    #  path.display
    # }) # renderText
   })# observeEvent r1_files

   observeEvent(input[[paste0('setup_grp',i,'_r2_files')]],{
    file = input[[paste0('setup_grp',i,'_r2_files')]]
    new.path=as.character(parseFilePaths(root = volumes,file)$datapath)
    new.path=gsub('^/filepath/','',new.path)
    new.path=gsub('^/root/','',new.path)
    text.path=input[[paste0('setup_grp',i,'_r2_filepaths')]]
    text.path=unlist(strsplit(text.path,"\n"))
    # combine new.path and paths in textarea
    setup.grp.r2.file.lst[[paste0('grp',i)]] <- sort(unique(c(new.path,text.path)))
    path.display=gsub('^/filepath/','',setup.grp.r2.file.lst[[paste0('grp',i)]])
    path.display=gsub('^/root/','',path.display)
    if(sum(path.display=="")>0){
     rem.id=which(path.display=="")
     path.display=path.display[-rem.id]
    }
    
    path.display <- paste0(path.display, collapse="\n")
    updateTextAreaInput(session, paste0('setup_grp',i,'_r2_filepaths'),value=path.display)
    # output[[paste0('setup.grp',i,'.r1.filepaths')]]<- renderText({
    #  path.display=gsub('^/filepath/','',setup.grp.r1.file.lst[[paste0('grp',i)]])
    #  path.display=gsub('^/root/','',path.display)
    #  path.display <- paste0(path.display[-1], "\n")
    #  path.display
    # }) # renderText
   })# observeEvent r1_files
    } #function
 )#lapply


    # # # get current content of textareainput
    # setup.grp.r1.file.lst[[paste0('grp',i)]] <- input[[paste0('setup_grp',i,'_r1_filepaths')]]
    # output[[paste0('setup.grp',i,'.r1.filepaths')]]<- renderText({
    #  path.display=gsub('^/filepath/','',setup.grp.r1.file.lst[[paste0('grp',i)]])
    #  path.display=gsub('^/root/','',path.display)
    #  path.display <- paste0(path.display[-1], "\n")
    #  path.display
    # }) # renderText
   
 # select project directory
 shinyDirChoose(input, "setup_proj_dir", root=volumes)

react.setup.proj.dir <- reactive({
  proj.dir = input$setup_proj_dir
  proj.path=parseDirPath(roots = volumes,proj.dir)
  proj.path=paste0(proj.path,"/")
 })

output$fileExists <- reactive({
 proj.dir <- react.setup.proj.dir ()
 file_exists=file.exists(file.path(proj.dir,"samples.txt"))
 file_exists
})
outputOptions(output, 'fileExists', suspendWhenHidden=FALSE)

 output$setup.proj.dir.filepaths <- renderText({
  projdir <- react.setup.proj.dir()
  path.display=gsub('^/filepath/','',projdir)
  path.display=gsub('^/root/','',path.display)
  path.display
 })
 
 observeEvent(input$setup.load.samples,{
  projdir <- react.setup.proj.dir()
  # load existing samples treating "NA" as regulating string
  df = read.table(file.path(projdir,"samples.txt"),na.strings=character(0),stringsAsFactors = FALSE)
  updateTextInput(session,inputId="setup.email",value=df[1,5])
  setup.btn <- setup.value()
  # number of rows to add
  add.row <- dim(df)[1]-setup.btn
  for (i in (1:dim(df)[1])){
   # ids for entries in each row
   setup.id.grp.name <- paste0("setup.grp",i,".name")
   setup.id.ctrl.name <- paste0("setup.grp",i,".ctrl.name")
   setup.id.rep.name <- paste0("setup.grp",i,".rep.name")
   setup.id.r1.files <- paste0("setup_grp",i,"_r1_files")
   setup.id.r2.files <- paste0("setup_grp",i,"_r2_files")
   setup.id.r1.filepath <- paste0("setup_grp",i,"_r1_filepaths")
   setup.id.r2.filepath <- paste0("setup_grp",i,"_r2_filepaths")
   # formatting path read from table
   # remove hostfilepath 
   if(hostfilepath!="/"){
    setup.df.path.display.r1=gsub(hostfilepath,'',df[i,6])
    setup.df.path.display.r1 <- gsub(",","\n",setup.df.path.display.r1)
    setup.df.path.display.r2=gsub(hostfilepath,'',df[i,7])
    setup.df.path.display.r2 <- gsub(",","\n",setup.df.path.display.r2)
   }
   if(hostfilepath=="/"){
    setup.df.path.display.r1=gsub("^/",'',df[i,6])
    setup.df.path.display.r1 <- gsub(",/","\n",setup.df.path.display.r1)
    setup.df.path.display.r2 <- gsub("^/","",df[i,7])
    setup.df.path.display.r2 <- gsub(",/","\n",setup.df.path.display.r2)
   }
   # add new rows if needed
  if(i > (dim(df)[1]-add.row)){
   # update the reactiveVal
   setup.value(i)
  # id for new div
  setup.div.id <- paste0("setup_div_",i)

  if(i %% 2 == 0) row.color="#FFFFFF"
  if(i %% 2 != 0) row.color="#f9f9f9"
  insertUI(
   selector = "#setup-placeholder",
   ui = tags$div(
    tags$style(HTML(paste0("#setup_grp",i,"_r1_filepaths { width: 300px; overflow-x: auto; white-space: pre;}"))),
    tags$style(HTML(paste0("#setup_grp",i,"_r2_filepaths { width: 300px; overflow-x: auto; white-space: pre;}"))),
    fluidRow(style = paste0("background-color:",row.color),
             column(2,
                    column(2,
                           tags$label(paste0(i,"."), style = "padding-top: 30px;"),
                           style = "padding-right: 0px;padding-left:0px"),
                    column(10,textAreaInput(inputId=setup.id.grp.name,
                                        label = "",value=df[i,1]),
                           style = "padding-right: 0px;padding-left:0px")
             ),
             column(2,textAreaInput(inputId=setup.id.ctrl.name,
                                label = "",value=df[i,2] )),
             column(2,textInput(inputId=setup.id.rep.name,
                                label = "",value=df[i,3] )),
             br(),
             column(3,shinyFilesButton(setup.id.r1.files,
                                       label = "Select R1 fastq files",
                                       title=paste("Please select Group",i,"R1 fastq files"),
                                       multiple=TRUE ),
                    textAreaInput(setup.id.r1.filepath,label="",width='150px',height='100px',
                                  value=setup.df.path.display.r1)),
             column(3,shinyFilesButton(setup.id.r2.files,
                                       label = "Select R2 fastq files",
                                       title=paste("Please select Group",i,"R1 fastq files"),
                                       multiple=TRUE),
                    textAreaInput(setup.id.r2.filepath,label="",width='150px',height='100px',
                                  value=setup.df.path.display.r2)),
             
    ),#fluidrow
    id = setup.div.id
   )
  )#insertUI
  setup.inserted.div <<- c(setup.inserted.div, setup.div.id)
  } # if(i > (dim(df)[1]-add.row))
   # fill out rows
   updateTextAreaInput(session,setup.id.grp.name,value = df[i,1])
   updateTextAreaInput(session,setup.id.ctrl.name,value = df[i,2])
   updateTextInput(session,setup.id.rep.name,value = df[i,3])
   updateTextAreaInput(session,setup.id.r1.filepath,value = setup.df.path.display.r1)
   updateTextAreaInput(session,setup.id.r2.filepath,value = setup.df.path.display.r2)
  }# for (i in (1:dim(df)[1]))
 })
 
 # select other genome files
 shinyFileChoose(input, "setup_genome_fa", root=volumes,
                 filetypes=c('', 'gz','fasta','fa'))
 shinyFileChoose(input, "setup_genome_gtf", root=volumes,
                 filetypes=c('', 'gz','GTF','gtf'))
 
 observeEvent(input$setup_genome_fa,{
  fa.file = input$setup_genome_fa
  fa.path=as.character(parseFilePaths(root = volumes,fa.file)$datapath)
  output$setup.genome.fa.filepaths <- renderText({
   path.display=gsub('^/filepath/','',fa.path)
   path.display=gsub('^/root/','',path.display)
   path.display
  })
 })
 
 observeEvent(input$setup_genome_gtf,{
  file = input$setup_genome_gtf
  path=as.character(parseFilePaths(root = volumes,file)$datapath)
  output$setup.genome.gtf.filepaths <- renderText({
   path.display=gsub('^/filepath/','',path)
   path.display=gsub('^/root/','',path.display)
   path.display
  })
 })
 
 output$err_msg_genome <- renderText({
  genome.fa = react.setup.genome.fa()
  genome.name=react.setup.genome.name()
  genome=react.setup.genome()
  genome.gtf=react.setup.genome.gtf()
  if(genome=="other"){
  shiny::validate(
   shiny::need(genome.name != '', "You must enter a genome name"),
   shiny::need(genome.fa != '', 'You must select genome fasta file'),
   shiny::need(genome.gtf != '', "You must select a GTF file")
  )
  }
 })
 
 output$err_msg <- renderText({
  projdir=react.setup.proj.dir()
  grpname=react.setup.grp.name()
  ctrlname=react.setup.grp.ctrl.name()
  repname=react.setup.grp.rep.name()
  nsamples=setup.value()
 shiny::validate(
  shiny::need(projdir != '', 'You must select a project folder'),
  shiny::need(!grepl("\\s",projdir), 'Project folder name cannot have spaces'),
  shiny::need(file.access(projdir,2) ==0, 'Please select a different project folder where you have write access'),
  shiny::need(length(ctrlname) == nsamples, 'You must enter control name'),
  shiny::need(length(repname)==nsamples,'You must enter replicate name'),
  shiny::need(nchar(repname)> 2,
              'You must use at least 3 characters for replicate name'),
  shiny::need(nchar(ctrlname)> 2,
              'You must use at least 3 characters for control name'),
  shiny::need(nchar(grpname)> 2,
              'You must use at least 3 characters for control name'),
  shiny::need(sum(grepl("na",ctrlname))==0,
              'You must use uppercase NA not lowercase na in ctrl name')
 )
 })
 
 observeEvent(input$setup.view.samples,{
  projdir=react.setup.proj.dir()
  req(paste0(projdir,"/samples.txt"))
  samplesContent <- readLines(paste0(projdir,"/samples.txt"))
  showModal(modalDialog(
   title="Samples.txt",
   verbatimTextOutput("samplesContent"),
   footer = actionButton("closeModal", "Close")
  ))
  output$samplesContent <- renderPrint({
   samplesContent
  })
  observeEvent(input$closeModal, {
    removeModal()
  })
 })
 
 observeEvent(input$setup.save.samples,{
  # create samples.txt
  projdir=react.setup.proj.dir()
  grpname=react.setup.grp.name()
  ctrlname=react.setup.grp.ctrl.name()
  repname=react.setup.grp.rep.name()
  spikename=rep("NA",length(grpname))
  email=react.setup.email()
  nsamples=setup.value()
  genome=react.setup.genome()
  genome.name = react.setup.genome.name()
  genome.fa=react.setup.genome.fa()
  genome.gtf=react.setup.genome.gtf()
  if(email==""){
   email=rep("NA",length(grpname))
  }else{
   email=rep(email,length(grpname))
  }
  print("projdir")
  print(projdir)
  shiny::req(projdir!='')
  shiny::req(length(grpname)==nsamples)
  shiny::req(length(ctrlname)==nsamples)
  shiny::req(length(repname)==nsamples)
  shiny::req(file.access(projdir,2) ==0)
  
  # getting the r1 and r2 fastq for all groups
  print("here")
  # function to remove empty paths
  rem_empty_path <- function(input.name) {
   filt.path = input[[input.name]]
   filt.path = unlist(strsplit(filt.path,"\n"))
   print("filt.path")
   print(filt.path)
   empty.id = which(filt.path=="")
   print("empty.id")
   print(empty.id)
   if(length(empty.id)>0) filt.path = filt.path[-empty.id]
   print("filt.path after rem empty")
   print(filt.path)
  }
  print("input[[setup_grp1_r1_filepaths)]]")
  print(input[[paste0("setup_grp",1,"_r1_filepaths")]])
  r1_fastq_lst=lapply(1:length(grpname),function(x) rem_empty_path(paste0("setup_grp",x,"_r1_filepaths")))
  r2_fastq_lst=lapply(1:length(grpname),function(x) rem_empty_path(paste0("setup_grp",x,"_r2_filepaths")))
  print("r1_fastq_lst")
  print(r1_fastq_lst)
  # sort and unique to make sure the correct pairing of r1 and r2
  r1_fastq_lst=lapply(1:length(grpname),
                      function(x) sort(unique(r1_fastq_lst[[x]])))
  r2_fastq_lst=lapply(1:length(grpname),
                      function(x) sort(unique(r2_fastq_lst[[x]])))
  # changing to hostpath
  r1_fastq=lapply(1:length(grpname),
                  function(x) paste(paste0(hostfilepath,r1_fastq_lst[[x]]),collapse=","))
  r2_fastq=lapply(1:length(grpname),
                  function(x) paste(paste0(hostfilepath,r2_fastq_lst[[x]]),collapse=","))
  print("r1_fastq")
  print(r1_fastq)
  df <- data.frame(
   grpname=grpname,
   ctrlname=ctrlname,
   repname=repname,
   spikename=spikename,
   email=email,
   path_to_r1_fastq=unlist(r1_fastq),
   path_to_r2_fastq=unlist(r2_fastq))
  print("file.path(projdir,samples.txt)")
  print(file.path(projdir,"samples.txt"))
  if(file.exists(file.path(projdir,"samples.txt"))){
   exist.samples.md5 = tools::md5sum(file.path(projdir,"samples.txt"))
  }
  write.table(paste0("#Groupname\tControlname\tReplicatename\tspikename\temail",
                     "\tpath_to_r1_fastq\tpath_to_r2_fastq"),
              file=file.path(projdir,"tmp_samples.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(df,file=file.path(projdir,"tmp_samples.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE,
              append = TRUE)
  tmp.samples.md5=tools::md5sum(file.path(projdir,"tmp_samples.txt"))
  # rename tmp.samples.txt to samples.txt if the content changes
  # this is to allow nextflow to resume
  if(file.exists(file.path(projdir,"samples.txt"))){
  if(!tmp.samples.md5 == exist.samples.md5)
   file.rename(file.path(projdir,"tmp_samples.txt"),file.path(projdir,"samples.txt"))
  }else{
   # there is no existing samples.txt just rename
   file.rename(file.path(projdir,"tmp_samples.txt"),file.path(projdir,"samples.txt"))
  }
 })
 
 observeEvent(input$setup.run.analysis,{
  # create samples.txt
   projdir=react.setup.proj.dir()
   grpname=react.setup.grp.name()
   ctrlname=react.setup.grp.ctrl.name()
   repname=react.setup.grp.rep.name()
   spikename=rep("NA",length(grpname))
   email=react.setup.email()
   nsamples=setup.value()
   genome=react.setup.genome()
   genome.name = react.setup.genome.name()
   genome.fa=react.setup.genome.fa()
   genome.gtf=react.setup.genome.gtf()
   if(email==""){
    email=rep("NA",length(grpname))
   }else{
    email=rep(email,length(grpname))
   }
   shiny::req(projdir!='')
   shiny::req(length(grpname)==nsamples)
   shiny::req(length(ctrlname)==nsamples)
   shiny::req(length(repname)==nsamples)
   shiny::req(file.access(projdir,2) ==0)
   
   # disable the run analysis button after it's clicked
   shinyjs::disable("setup.run.analysis")
   shinyjs::show("text1")
   
   # getting the r1 and r2 fastq for all groups
   # function to remove empty paths
   rem_empty_path <- function(input.name) {
    filt.path = input[[input.name]]
    filt.path = unlist(strsplit(filt.path,"\n"))
    print("filt.path")
    print(filt.path)
    empty.id = which(filt.path=="")
    print("empty.id")
    print(empty.id)
    if(length(empty.id)>0) filt.path = filt.path[-empty.id]
    print("filt.path after rem empty")
    print(filt.path)
	}
   print("input[[setup_grp1_r1_filepaths)]]")
   print(input[[paste0("setup_grp",1,"_r1_filepaths")]])
   r1_fastq_lst=lapply(1:length(grpname),function(x) rem_empty_path(paste0("setup_grp",x,"_r1_filepaths")))
   r2_fastq_lst=lapply(1:length(grpname),function(x) rem_empty_path(paste0("setup_grp",x,"_r2_filepaths")))
   print("r1_fastq_lst")
   print(r1_fastq_lst)
   # sort and unique to make sure the correct pairing of r1 and r2
   r1_fastq_lst=lapply(1:length(grpname),
          function(x) sort(unique(r1_fastq_lst[[x]])))
   r2_fastq_lst=lapply(1:length(grpname),
           function(x) sort(unique(r2_fastq_lst[[x]])))
   # changing to hostpath
   r1_fastq=lapply(1:length(grpname),
                   function(x) paste(paste0(hostfilepath,r1_fastq_lst[[x]]),collapse=","))
   r2_fastq=lapply(1:length(grpname),
                   function(x) paste(paste0(hostfilepath,r2_fastq_lst[[x]]),collapse=","))
   print("r1_fastq")
   print(r1_fastq)
   df <- data.frame(
    grpname=grpname,
    ctrlname=ctrlname,
    repname=repname,
    spikename=spikename,
    email=email,
    path_to_r1_fastq=unlist(r1_fastq),
    path_to_r2_fastq=unlist(r2_fastq))
   
   # only create samples.txt if the content changes
   # this is to allow nextflow to resume
      if(file.exists(file.path(projdir,"samples.txt"))){
    exist.samples.md5 = tools::md5sum(file.path(projdir,"samples.txt"))
   }
   write.table(paste0("#Groupname\tControlname\tReplicatename\tspikename\temail",
                      "\tpath_to_r1_fastq\tpath_to_r2_fastq"),
               file=file.path(projdir,"tmp_samples.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
   write.table(df,file=file.path(projdir,"tmp_samples.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE,
               append = TRUE)
   tmp.samples.md5=tools::md5sum(file.path(projdir,"tmp_samples.txt"))

   if(file.exists(file.path(projdir,"samples.txt"))){
   if(!tmp.samples.md5 == exist.samples.md5)
    file.rename(file.path(projdir,"tmp_samples.txt"),file.path(projdir,"samples.txt"))
   }else{
    # just rename there is no existing samples.txt
    file.rename(file.path(projdir,"tmp_samples.txt"),file.path(projdir,"samples.txt"))
   }
  #system("echo 'sbatch --help' > /hostpipe")
  # getting genome options
  if(genome=="other"){
   shiny::req(genome.fa != "")
   shiny::req(genome.gtf!="")
   shiny::req(genome!="")

   # convert path to hostpath
   if(file.exists("/filepath")){
    hostfa=gsub('/filepath/',hostfilepath,genome.fa)
    hostgtf=gsub('/filepath/',hostfilepath,genome.gtf)
   }else{ 
    hostfa=gsub('/root/','/',genome.fa)
    hostgtf=gsub('/root/','/',genome.gtf)
   }
   options=paste0("genome=",genome.name," ref_fa=",hostfa," ref_gtf=",hostgtf)
   if(run_nextflow)
	  options=paste0("--genome=",genome.name," --fasta=",hostfa," --gtf=",hostgtf)
  }else{
   options=paste0("genome=",genome)
   if(run_nextflow){
    options=paste0("--genome=",genome)
    if(genome=="hg19"){
	    # search for fasta file in /ref/hg19
     # /ref is bound by singularity to img.dir
     fa.name=list.files(path="/ref/hg19",pattern="\\.(fasta|fa)$")
	    gtf.name=list.files(path="/ref/hg19",pattern="\\.(gtf)$")
	    print("fa.name")
	    print(fa.name)
	    # options given to nextflow need to have full path including img.dir
	    options=paste0(options," --fasta=",img.dir,"/ref/hg19/",fa.name[1],
	                   " --gtf=",img.dir,"/ref/hg19/",gtf.name[1])
    }else if(genome=="hg38"){
     # search for fasta file in /ref/hg38
     # /ref is bound by singularity to img.dir
     fa.name=list.files(path="/ref/hg38",pattern="\\.(fasta|fa)$")
     gtf.name=list.files(path="/ref/hg38",pattern="\\.(gtf)$")
     print("fa.name")
     print(fa.name)
     # options given to nextflow need to have full path including img.dir
     options=paste0(options," --fasta=",img.dir,"/ref/hg38/",fa.name[1],
                    " --gtf=",img.dir,"/ref/hg38/",gtf.name[1])
    }
    }
   }
   print("options") 
   print(options)
  if(input$setup.batch.adjust==TRUE){ batch_adjust="yes"
  }else{ batch_adjust="no"}
   
   if(input$setup.debug==TRUE){ run="debug"
   }else{ run=""}
  if(!run_nextflow){
   options=paste0(options," time=",input$setup.time," batch_adjust=",batch_adjust,
                 " ncpus_trim=",input$setup.ncpus.trim," ncpus_star=",input$setup.ncpus.star,
                 " run=",run)
  }else{
	  options=paste0(options," --batch_adjust=",batch_adjust)
  }
  # converting project dir to host path
  if(file.exists("/filepath")){
   hostprojdir=gsub('/filepath/',hostfilepath,projdir)
  }else{ 
   hostprojdir=gsub('/root/','/',projdir)
  }
   if(!run_nextflow){
   system(paste0("echo 'cd ",hostprojdir,"&& bash ",img.dir,"/scripts/run_rnaseq_full.sh ",options,
                " &> run_rnaseq_full.out' > /hostpipe"))
   }else{
    # create nextflow config
    # read variables defined in nextflow_config_var.config
    text.nf=readLines("/scripts/nextflow_config_var.config")
    for (i in seq_along(text.nf)){
     if(grepl("#",text.nf[i]))
      assign(strsplit(text.nf[i+1],"=")[[1]][1],strsplit(text.nf[i+1],"=")[[1]][2])
    }
    # convert time to nextflow format
    nf.setup.time=""
    tmp.setup.time=unlist(strsplit(input$setup.time,"-"))
    if(length(tmp.setup.time) > 1){
     nf.setup.time=paste0(tmp.setup.time[1],"d")
     tmp.setup.time=tmp.setup.time[2]
    }
     tmp.setup.time=unlist(strsplit(tmp.setup.time,":"))
     # pad zeros if not defined
     while(length(tmp.setup.time)<3){
      tmp.setup.time[length(tmp.setup.time)+1]=0
     }
    nf.setup.time=paste(nf.setup.time,
                         paste0(tmp.setup.time,c('h','m', 's'),collapse=" "))
     nf.setup.time=trimws(nf.setup.time)
     print("email")
     print(email)
     # remove extra quote
     addtl_opt <- gsub("\"","",addtl_opt)
     print(addtl_opt)
     if(email[1]=="NA")
      email=""
    nf.config=paste0(
    "params {
      help= false
      version = false
      monochrome_logs= false
    }\n\n",
    "process {\n",
      "// container path defined as absolute path to <img.dir> rnaseq-pipe-container.sif\n",
      " container = 'file:///",img.dir,"/rnaseq-pipe-container.sif'\n",
      " time = '",nf.setup.time,"'\n",
      " queue = '",general_partition,"'\n",
      " cpus = 1\n",
      " email = '",email,"'\n",
      " clusterOptions = '",addtl_opt,"'\n\n",

      "withLabel: hi_mem {\n",
      "  queue= '",high_mem_partition,"'\n",
      "  memory= '",high_mem,"g'\n",
      "}\n\n",

      "withLabel: TRIM_FASTQC {\n",
      "  cpus = ",input$setup.ncpus.trim,"\n",
      "}\n\n",
      
      "withLabel: short_time {\n",
      "  time = '30min'\n",
      "}\n\n",
      
      "withLabel: hi_cpus {\n",
      "  cpus = 15\n",
      "}\n\n",
      
      "// used for star index\n",
      "withLabel: star {\n",
      "   queue = '",high_mem_partition,"'\n",
      "   cpus = ",input$setup.ncpus.star,"\n",
      "}\n\n",
      "withLabel: very_hi_mem {\n",
      "  queue = '",high_mem_partition,"'\n",  
      "  memory = '",very_high_mem,"g'\n",
      "}\n\n",
      "// used for star pass1 and pass2\n",
      "withLabel: hi_mem_cpus {\n",
      " queue = '",high_mem_partition,"'\n",
      " cpus = ",input$setup.ncpus.star,"\n",
    "}\n\n",
    
    "}\n",
    "singularity.enabled = true\n",
    "singularity.autoMounts = true\n\n",
    
    "process.executor = '",executor,"'\n",
    "// check job status in all partition in case job get reassigned to diff partition\n",
    "executor.queueGlobalStatus = true\n",
    "dag.overwrite = true\n",
    "report.overwrite = true\n"
    )
    # get mount point to bind to singularity
    mount_point=strsplit(hostprojdir,"/")[[1]][2];
    mount_point_str=paste0(",/",mount_point,":/",mount_point,"\"")

    if(!genome=="other"){
     nf.config=paste0(nf.config,"singularity.runOptions = \"--bind ",
                      img.dir,"/scripts:/scripts,",
     #                 hostprojdir,":/mnt,",img.dir,"/ref:/ref",",/gpfs0:/gpfs0\"")
     #                 hostprojdir,":/mnt,",img.dir,"/ref:/ref\"")
			hostprojdir,":/mnt,",img.dir,"/ref:/ref",mount_point_str)

    }else{
     hostfolderfa=dirname(hostfa)
     nf.config=paste0(nf.config,'singularity.runOptions = \"--bind ',
                      img.dir,'/scripts:/scripts,',hostfolderfa,":/ref,",
     #                 hostprojdir,':/mnt,/gpfs0:/gpfs0\"')
     #                 hostprojdir,':/mnt\"')
		       hostprojdir,":/mnt,",mount_point_str)
    }
    print("nf.config")
    print(nf.config)
    writeLines(nf.config,"/img_dir/nextflow.config")
    # cmd=paste0("echo 'echo \"",nf.config," >> ",img.dir,"/nextflow.config",
    #            " && ml load nextflow/22.10.6.5843 && nextflow run -resume ",img.dir,
    #            "/main.nf -ansi-log false --inputdir=",hostprojdir," ",options,
    #            " &> run_rnaseq_full.out' > /hostpipe")
    if(is.na(load_nextflow))
     load_nextflow_cmd=""
    else
     load_nextflow_cmd=paste0(load_nextflow," && ")
    cmd=paste0("echo '",load_nextflow_cmd,"cd ", hostprojdir, 
               " && nextflow run -resume ",img.dir,
               "/main.nf -ansi-log true -with-report ",
               hostprojdir,"run_rnaseq_full.html --inputdir=",
               hostprojdir," ",options," &> ",hostprojdir,"run_rnaseq_full.out' > /hostpipe")
    print(cmd)
    system(cmd)
   }
})
 
 #### log tab #####
 # update log file list
 observeEvent(input$logtab.refresh.log.path,{
  # list of files in log directory from most recent
  projdir <- react.setup.proj.dir()
  # copy log (i.e. *.out) files to log folder to view
  files.to.copy <- list.files(projdir,pattern="\\.out$",full.names=TRUE)
  # add run_shiny_analysis.out
  if(file.exists("run_shiny_analysis.out")) 
   files.to.copy <- c(files.to.copy,"run_shiny_analysis.out")
  if(run_nextflow)
   files.to.copy <- c(files.to.copy,file.path(projdir,".nextflow.log"))
  print("files.to.copy")
  print(files.to.copy)
  print("projdir/outputs/logs")
  print(paste0(projdir,"outputs/logs"))
  # Check if the logs directory exists; if not, create it
  dir_path=paste0(projdir,"outputs/logs")
  if (!dir.exists(dir_path)) {
   dir.create(dir_path, recursive = TRUE)
   cat("Directory created:", dir_path, "\n")
  } else {
   cat("Directory already exists:", dir_path, "\n")
  }
  #file.copy.msg=file.copy(files.to.copy,paste0(projdir,"outputs/logs"),overwrite=TRUE)
  # modified to use sapply to avoid error when copying multiple files
  file.copy.msg=sapply(files.to.copy, 
                       function(x) file.copy(x,paste0(projdir,"outputs/logs"),overwrite = TRUE))
  files=list.files(paste0(projdir,"outputs/logs"),
                   pattern = "\\.txt$|\\.out$|\\.log$",full.names = TRUE,all.files = TRUE)
  print("files")
  print(files)
  # Sort files, placing filenames starting with "run_" at the top
  sorted_files <- c(sort(files[startsWith(files, "run_")]), 
                    sort(files[!startsWith(files, "run_")]))
  log.lst = basename(sorted_files)
  if(file.exists(paste0(projdir,"outputs/logs/run_rnaseq_full.out"))){
   updateSelectInput(session,"logtab.log.path",choices=log.lst,selected = "run_rnaseq_full.out")
  }else{
   updateSelectInput(session,"logtab.log.path",choices=log.lst,selected="")
  }
     # display file content in UI
   output$logtab_log_content <- renderPrint({ 
    shiny::req(input$logtab.log.path!="")
    cat(react.logtab_log_content(),sep="\n")
   })
  # list only failed *.out file
  files=list.files(paste0(projdir,"outputs/logs"),pattern = "\\.out$",full.names = TRUE)
  # Filter files containing the word "FAILED"
  failed.files <- lapply(files, function(file) {
   if (any(grepl("FAILED", tail(readLines(file),2), ignore.case = FALSE))) {
    return(file)
   }
  })
  # Filter out NULL values (files without the word "FAILED")
  failed.files <- unlist(Filter(Negate(is.null), failed.files))
  if(is.null(failed.files)){
   log.fail.lst=""
  }else{
   # get the first FAILED
   file_info <- file.info(failed.files)
   sorted_files <- failed.files[order(as.numeric(file_info$mtime),decreasing=FALSE)]
   log.fail.lst = basename(sorted_files)
  }
  updateSelectInput(session,"logtab.fail.log.path",choices=log.fail.lst)
 
  })#observEvent refresh
 
 react.logtab_log_content <- reactive({
    projdir <- react.setup.proj.dir()
    path <- paste0(projdir,"outputs/logs/",input$logtab.log.path)
  content <- readLines(path)
  return(content)
 })
 
 react.logtab.fail.log.content <- reactive({
  projdir <- react.setup.proj.dir()
  if(input$logtab.fail.log.path=="NA") return("Nothing failed!")
  path <- paste0(projdir,"outputs/logs",input$logtab.fail.log.path)
  content <- readLines(path)
  return(content)
 })
 
 # display file content in UI
 output$logtab.fail.log.content <- renderPrint({
  shiny::req(input$logtab.fail.log.path!="")
  react.logtab.fail.log.content()
 })
 
 # Observe the global variable and add/remove a tab accordingly
 observe({
  if (run_nextflow) {
   # Add a new tab if the global variable is TRUE
   insertTab(inputId = "tabset",
             tabPanel(id="nfreporttab", "Nextflow report", fluid=TRUE,
                      htmlOutput("nfreport")
             ),
             target = "Log",  # Insert the tab after "Log Tab"
             position = "after")
   # hide run debug
   hide("setup.debug")
  }
  
 })
 
 #### Nextflow report tab #####
 output$nfreport <- renderUI({
  projdir <- react.setup.proj.dir()
  # Check if the file exists
  report_path <- file.path(projdir,"run_rnaseq_full.html")
  if (!file.exists(report_path)) {
    return(
      tags$div(
        style = "color: blue; font-weight: bold; padding: 20px;",
    tags$p("Nextflow report not found. Possible reasons:"),
    tags$ul( style = "color: black; font-weight:normal; padding: 20px;",
      tags$li("The pipeline may not have finished."),
      tags$li("The project folder might be incorrect."),
      tags$li("The analysis may not have been run using RNA-SeqEZPZ-NF.")
    )
      )
    )
  }

  # path for nfreport html to be referred to in iframe
  # using unique alias so the path will be updated when 
  # projdir is updated
  uid <- paste0("nfreport_",as.integer(Sys.time()))
  addResourcePath(uid,projdir)
  tags$iframe(seamless="seamless",
              #src="nfreport/run_rnaseq_full.html",
              src=file.path(uid,"run_rnaseq_full.html"),
              width="100%",
              height="1000")
 })
 
 #### QC tab #####
 output$multiqc <- renderUI({
  projdir <- react.setup.proj.dir()
 # Check if the file exists
  report_path <- file.path(projdir,"outputs","fastqc_rslt","multiqc_report.html")
  if (!file.exists(report_path)) {
    return(
      tags$div(
        style = "color: blue; font-weight: bold; padding: 20px;",
    tags$p("QC report not found. Possible reasons:"),
    tags$ul( style = "color: black; font-weight:normal;	padding: 20px;",
      tags$li("The pipeline may not have finished."),
      tags$li("The project folder might be incorrect.")
    )	    
      )
    )
  }
  # path for fastqc report html to be referred to in iframe
  # using unique alias so the path will be updated when 
  # projdir is updated
  uid <- paste0("fastqc_rslt_",as.integer(Sys.time()))
  addResourcePath(uid,file.path(projdir,"outputs","fastqc_rslt"))
  #addResourcePath("fastqc_rslt", paste0(projdir,"outputs/fastqc_rslt"))
  tags$iframe(seamless="seamless",
              #src="fastqc_rslt/multiqc_report.html",
              src=file.path(uid,"multiqc_report.html"),
              width="100%",
              height="1000")
 })
 
 ## Differential genes analysis report ###
 output$diff_report <- renderUI({
  projdir <- react.setup.proj.dir()
  # Check if the file exists
  report_path <- file.path(projdir,"outputs","diff_analysis_rslt","RNA-seq_differential_analysis_report.html")
  if (!file.exists(report_path)) {
    return(
      tags$div(
        style = "color: blue; font-weight: bold; padding: 20px;",
    tags$p("Differential genes analysis not found. Possible reasons:"),
    tags$ul( style = "color: black; font-weight:normal;	padding: 20px;",
      tags$li("The pipeline may not have finished."),
      tags$li("The project folder might be incorrect.")
    )	    
      )
    )
  }
  # path for differential report html to be referred to in iframe
  # using unique alias so the path will be updated when 
  # projdir is updated
  uid <- paste0("diff_report_path_",as.integer(Sys.time()))
  #addResourcePath("diff_report_path",paste0(projdir,"outputs/diff_analysis_rslt/"))
  addResourcePath(uid,file.path(projdir,"outputs","diff_analysis_rslt"))
  print("uid")
  print(uid)
  tags$iframe(seamless="seamless", 
              #src= "diff_report_path/RNA-seq_differential_analysis_report.html",
              src= file.path(uid,"RNA-seq_differential_analysis_report.html"),
              width="100%", 
              height=800)
 })
 
 #### processing plots tab #####
 
 react.tab3.data.loc <- reactive({
  projdir <- react.setup.proj.dir()
  data.loc = paste0(projdir,"outputs/")
 })
 
 react.tab3.out.loc <- reactive({
   data.loc <- react.tab3.data.loc()
   shiny::req(dir.exists(data.loc))
  # create directory for shiny outputs
  dir.create(paste0(data.loc, "shiny_outputs"))
  out.loc=paste0(data.loc,"shiny_outputs/")
 })
 
 # initialize values for tab
react.tab3.rdata <- reactive({
 req(input$tabset=="Plots")
 projdir <- react.setup.proj.dir()
 data.loc <- react.tab3.data.loc()
 rdata.loc1=paste0(data.loc,"diff_analysis_rslt/RNA-seq_differential_analysis.RData")
 rdata.loc2=paste0(data.loc,"diff_analysis_rslt/RNA-seq differential analysis.RData")
 shiny::req(file.exists(rdata.loc1)||file.exists(rdata.loc2))
 # filename was changed. Make sure can read both
 if(file.exists(rdata.loc1)){
  load(rdata.loc1)
   } else if (file.exists(rdata.loc2)){
  load(rdata.loc2)
   }
 return(out.DESeq2)
})

 react.tab3.groups.lst <- reactive({
  out.DESeq2 <- react.tab3.rdata()
  groups= data.frame(
  var = names(out.DESeq2$results),
  names = names(out.DESeq2$results)
 )
 # List of choices for groupnames
 groups.lst <- as.list(groups$names)
 # Name it
 names(groups.lst) <- groups$var
 groups.lst
})

 observe({
  groups.lst <- react.tab3.groups.lst()
  updateSelectInput(session, "tab3.grp1.name",choices=groups.lst)
  updateTextAreaInput(session, "tab3.grp1.plot.title",value=groups.lst[[1]])
 })
 
 observe({
  grp.name <- react.tab3.grp.name()
  for (i in 1:length(grp.name)){
   updateTextAreaInput(session, paste0("tab3.grp",i,".plot.title"),value=grp.name[i])
   }# for
 })#observe
 
 # vals will contain all plots and table grobs
 vals.plot <- reactiveValues(venn.up1=NULL,venn.up2=NULL,venn.dwn1=NULL,
                             venn.dwn2=NULL,hm.up=NULL,hm.dwn=NULL,
                             upset.up=NULL,upset.dwn=NULL,msig.mf=NULL,
                             msig.bp=NULL,msig.cc=NULL,msig.curate=NULL)
 vals.plot.volcano <- reactiveValues()
 
 # adding groups
 value <- reactiveVal(1)
 inserted <- c()		
 inserted.col <- c()
 inserted.fdr <- c()
 inserted.fc <- c()
 inserted.meanDiff <- c()
 # reactive value to plot gene enrichment
 react.val.gen.go <- reactiveVal(0)
 # adjust reactive value as gen.go button is pressed
 observeEvent(input$gen.go, {
  gen.go <- react.val.gen.go()
  gen.go <- react.val.gen.go(1)
 })
 
 observeEvent(input$insert_set, {
  groups.lst <- react.tab3.groups.lst()
  grpname <-react.tab3.grp.name()
  # remove selected choices (grpname) from list of choices (groups.lst)
  groups.lst <- setdiff(groups.lst,grpname)
  btn <- value() +1
  value(btn)
  if(btn %% 2 == 0) row.color="#FFFFFF"
  if(btn %% 2 != 0) row.color="#f9f9f9"
  id <- paste0("txt1_", btn)
  insertUI(
   selector = "#placeholder",
   ui = tags$div(
    fluidRow(style = paste0("background-color:",row.color),
     column(2,
    selectInput(inputId=paste0("tab3.grp",btn,".name"),
                label=(paste0("Group ",btn)), 
                selected=groups.lst[btn],
                choices=groups.lst),
    tags$style(    type = 'text/css',
                   ".selectize-input { word-wrap : break-word;}
                         .selectize-dropdown {word-wrap : break-word;} "
    )),
    column(2,textAreaInput(inputId=paste0("tab3.grp",btn,".plot.title"),
              label=(paste0("Group ",btn," label")), 
              value = groups.lst[btn])),
    column(2,textInput(inputId=paste0("tab3.color.grp",btn), 
              label = (paste0("Color for group ",btn)), 
              value = tab3.colors[btn])),
    column(2,numericInput(inputId=paste0("tab3.fdr.grp",btn), 
                          label = (paste0("FDR cut-off for group ",btn)), 
                          value = 0.05,step=0.01,min=0,max=1)),
    column(2,numericInput(inputId=paste0("tab3.fc.grp",btn), 
                 label = (paste0("Fold-change for group ",btn)), 
                 value = 1,min=1,step=0.5)),
    column(2,numericInput(inputId=paste0("tab3.meanDiff.grp",btn), 
                 label = (paste0("Difference for group ",btn)), 
                 value = 0,min=0,step=1))
    ),
    id = id
   )
  )
  inserted <<- c(inserted, id)
  
 })
 
 observeEvent(input$remove_set, {
  if(value()>1){
   btn <- value() -1
   value(btn)
   removeUI(
    selector = paste0("#",inserted[length(inserted)])
   )

   updateSelectInput(
    session,
    paste0("tab3.grp",length(inserted)+1 ,".name"),
    NULL,
    NA
   )
   updateTextAreaInput(
    session,
    paste0("tab3.grp", length(inserted)+1,".plot.title"),
    NULL,
    NA
   )

   updateTextInput(
    session,
    paste0("tab3.color.grp",length(inserted)+1),
    NULL,
    NA
   )
   updateNumericInput(
    session,
    paste0("tab3.fdr.grp",length(inserted)+1),
    NULL,
    NA
   )
   updateNumericInput(
    session,
    paste0("tab3.fc.grp",length(inserted)+1),
    NULL,
    NA
   )
   updateNumericInput(
    session,
    paste0("tab3.meanDiff.grp",length(inserted)+1),
    NULL,
    NA
   )
   
   inserted <<- inserted[-length(inserted)]
  }else{ btn <- 1
  value(btn)}
 })
 
 #### processing three groups #######
 react.grp.name <- reactive({
  c(input$grp1.name,input$grp2.name,input$grp3.name)
 })
 react.grp.plot.title <- reactive({
  c(input$grp1.plot.title,input$grp2.plot.title,input$grp3.plot.title)
 })
 
 react.fc.cutoff <- reactive({
  c(input$fc.grp1,input$fc.grp2,input$fc.grp3)
 })
 
 react.venn.cex <- reactive({
  input$venn.cex
 })	
 
 react.venn.pad <- reactive({
  input$venn.pad
 })		
 
 react.fdr.cutoff <- reactive({
  c(input$fdr.grp1,input$fdr.grp2,input$fdr.grp3)
 })
 
 react.venn.opts <- reactive({input$venn.opts})
 
 react.result <- reactive({
  grp.name=react.grp.name()
  grp.plot.title=react.grp.plot.title()
  idx=match(grp.name,names(out.DESeq2$results))
  results=out.DESeq2$results[idx]
 })
 
 gs.RNASeq <- reactive({
  grp.name=react.grp.name()
  grp.plot.title=react.grp.plot.title()
  results=react.result()
  gs.RNASeq <- max(sapply(results,function(x) sum(!is.na(x$padj))))
 })
 
 grp.up <- reactive ({
  grp.name=react.grp.name()
  grp.plot.title=react.grp.plot.title()
  results=react.result()
  fdr.co <- react.fdr.cutoff()
  fc.cutoff <- react.fc.cutoff()
  # up-regulated genes
  grp.up=list()
  for(i in 1:ngroup){
   grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
                                       results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) & 
                                       !is.na(results[[grp.name[i]]]$padj),]
  }
  grp.up
 })
 
 genes.lists.up <- reactive ({
  grp.up <- grp.up()
  grp.plot.title=react.grp.plot.title()
  genes.lists.up=lapply(grp.up,function(x) rownames(x))
  names(genes.lists.up)<-grp.plot.title
  genes.lists.up
 })
 
 gom.obj.up <- reactive({
  genes.lists.up <- genes.lists.up()
  gs.RNASeq <- gs.RNASeq()
  newGOM(genes.lists.up,genes.lists.up,genome.size=gs.RNASeq)
 })
 
 s4.up <- reactive({
  # plot euler venn
  s4<-genes.lists.up()
  names(s4)<-react.grp.plot.title()
  s4
 })
 
 grp.dwn <- reactive({
  grp.name=react.grp.name()
  grp.plot.title=react.grp.plot.title()
  results=react.result()
  fdr.co <- react.fdr.cutoff()
  fc.cutoff <- react.fc.cutoff()
  
  # down-regulated genes
  grp.dwn=list()
  for(i in 1:ngroup){
   grp.dwn[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
                                        results[[grp.name[i]]]$log2FoldChange<(-log2(fc.cutoff[i])) & 
                                        !is.na(results[[grp.name[i]]]$padj),]
  }
  grp.dwn
 })
 
 genes.lists.dwn <- reactive({
  grp.dwn <- grp.dwn()
  grp.plot.title=react.grp.plot.title()
  genes.lists.dwn <-
   lapply(grp.dwn,function(x) rownames(x))
  names(genes.lists.dwn)<-grp.plot.title
  genes.lists.dwn
 })
 
 gom.obj.dwn <- reactive({
  gs.RNASeq<-gs.RNASeq()
  genes.lists.dwn<-genes.lists.dwn()
  names(genes.lists.dwn)<-react.grp.plot.title()
  gom.obj.dwn <- newGOM(genes.lists.dwn,genes.lists.dwn,genome.size=gs.RNASeq)
  gom.obj.dwn
 })
 
 genes.lists <- reactive ({
  grp.name=react.grp.name()
  grp.plot.title=react.grp.plot.title()
  results=react.result()
  fdr.co <- react.fdr.cutoff()
  fc.cutoff <- react.fc.cutoff()
  grp.dwn <- grp.dwn ()
  grp.up <- grp.up ()
  # all regulated-genes
  grp=list()
  for(i in 1:ngroup){
   grp[[i]]=rbind(grp.up[[i]],grp.dwn[[i]])
  }
  
  genes.lists=lapply(grp,function(x) rownames(x))
  names(genes.lists)<-grp.plot.title
  genes.lists
 })
 
 gom.obj <- reactive({
  genes.lists <- genes.lists()
  gs.RNASeq <- gs.RNASeq()
  newGOM(genes.lists,genes.lists,genome.size=gs.RNASeq)
 })
 
 s4 <- reactive({
  # plot euler venn
  s4<-genes.lists()
  names(s4)<-react.grp.plot.title()
  s4
 })
 
 react.colors <- reactive({
  colors[1]=input$color.grp1
  colors[2]=input$color.grp2
  colors[3]=input$color.grp3
  colors
 })
 s4.dwn <- reactive({
  # plot euler venn
  s4<-genes.lists.dwn()
  names(s4)<-react.grp.plot.title()
  s4
 })
 
 data.up <- reactive({
  gom.obj.up <- gom.obj.up()
  grp.plot.title <- react.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj.up,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[c(2:3,6),]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj.up,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[c(2,3,6),]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 data.dwn <- reactive({
  gom.obj.dwn <- gom.obj.dwn()
  grp.plot.title <- react.grp.plot.title()
  ## heatmap down regulation
  data <- getMatrix(gom.obj.dwn,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[c(2,3,6),]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj.dwn,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[c(2,3,6),]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 data <- reactive({
  gom.obj <- gom.obj()
  grp.plot.title <- react.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[c(2:3,6),]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[c(2,3,6),]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 output$plot1 <- renderPlot({
  venn.opts=react.venn.opts()
  colors=react.colors()
  cex = react.venn.cex()
  pad = react.venn.pad()
  set.seed(1)
  s4 <- s4.up()
  venn.up <- plot_euler(s4,colors,cex,venn.opts,title="Up-regulation")
  s4 <- s4.dwn()
  venn.dwn <- plot_euler(s4,colors,cex,venn.opts,title="Down-regulation")
  grid.arrange(venn.up,venn.dwn,ncol=2,padding=unit(pad,"line"),
               bottom="",right="",left="",top="")
 })
 
 output$plot2 <- renderPlot({
   pad = react.venn.pad()
  data=data.up()
  hm.up <- heatmap_sigf_overlap(data)
  data=data.dwn()
  hm.dwn <- heatmap_sigf_overlap(data)
  grid.arrange(hm.up,hm.dwn,ncol=2,padding=unit(pad,"line"),
               bottom="",right="",left="",top="")
 })
 
 output$plot3 <- renderPlot({	
  venn.opts=react.venn.opts()
  colors=react.colors()
  cex = react.venn.cex()
  set.seed(1)
  s4 <- s4()
  plot_euler(s4,colors,cex,venn.opts,title="All regulated genes")
 })
 
 output$plot4 <- renderPlot({
  data=data()
  heatmap_sigf_overlap(data)
 })
 
 ######## processing two groups ############
 
 react.tab0.grp.name <- reactive({
  c(input$tab0.grp1.name,input$tab0.grp2.name)
 })
 react.tab0.grp.plot.title <- reactive({
  c(input$tab0.grp1.plot.title,input$tab0.grp2.plot.title)
 })
 
 react.tab0.venn.cex	<- reactive({
  input$tab0.venn.cex
 })
 
 react.tab0.venn.pad	<- reactive({
  input$tab0.venn.pad
 })
 
 react.tab0.fc.cutoff <- reactive({
  c(input$tab0.fc.grp1,input$tab0.fc.grp2)
 })
 
 react.tab0.fdr.cutoff <- reactive({
  c(input$tab0.fdr.grp1,input$tab0.fdr.grp2)
 })
 
 react.tab0.venn.opts <- reactive({input$tab0.venn.opts})
 
 react.tab0.result <- reactive({
  tab0.grp.name=react.tab0.grp.name()
  tab0.grp.plot.title=react.tab0.grp.plot.title()
  idx=match(tab0.grp.name,names(out.DESeq2$results))
  tab0.results=out.DESeq2$results[idx]
 })
 
 tab0.gs.RNASeq <- reactive({
  tab0.results=react.tab0.result()
  gs.RNASeq <- max(sapply(tab0.results,function(x) sum(!is.na(x$padj))))
 })
 
 tab0.grp.up <- reactive ({
  grp.name=react.tab0.grp.name()
  grp.plot.title=react.tab0.grp.plot.title()
  results=react.tab0.result()
  fdr.co <- react.tab0.fdr.cutoff()
  fc.cutoff <- react.tab0.fc.cutoff()
  # up-regulated genes
  grp.up=list()
  for(i in 1:tab0.ngroup){
   grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
                                       results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) & 
                                       !is.na(results[[grp.name[i]]]$padj),]
  }
  grp.up
 })
 
 tab0.genes.lists.up <- reactive ({
  grp.up <- tab0.grp.up()
  grp.plot.title=react.tab0.grp.plot.title()
  genes.lists.up=lapply(grp.up,function(x) rownames(x))
  names(genes.lists.up)<-grp.plot.title
  genes.lists.up
 })
 
 tab0.gom.obj.up <- reactive({
  genes.lists.up <- tab0.genes.lists.up()
  gs.RNASeq <- tab0.gs.RNASeq()
  newGOM(genes.lists.up,genes.lists.up,genome.size=gs.RNASeq)
 })
 
 tab0.s4.up <- reactive({
  # plot euler venn
  s4<-tab0.genes.lists.up()
  names(s4)<-react.tab0.grp.plot.title()
  s4
 })
 
 tab0.grp.dwn <- reactive({
  grp.name=react.tab0.grp.name()
  grp.plot.title=react.tab0.grp.plot.title()
  results=react.tab0.result()
  fdr.co <- react.tab0.fdr.cutoff()
  fc.cutoff <- react.tab0.fc.cutoff()
  
  # down-regulated genes
  grp.dwn=list()
  for(i in 1:tab0.ngroup){
   grp.dwn[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
                                        results[[grp.name[i]]]$log2FoldChange<(-log2(fc.cutoff[i])) & 
                                        !is.na(results[[grp.name[i]]]$padj),]
  }
  grp.dwn
 })
 
 tab0.genes.lists.dwn <- reactive({
  grp.dwn <- tab0.grp.dwn()
  grp.plot.title=react.tab0.grp.plot.title()
  genes.lists.dwn <-
   lapply(grp.dwn,function(x) rownames(x))
  names(genes.lists.dwn)<-grp.plot.title
  genes.lists.dwn
 })
 
 tab0.gom.obj.dwn <- reactive({
  gs.RNASeq<-tab0.gs.RNASeq()
  genes.lists.dwn<-tab0.genes.lists.dwn()
  names(genes.lists.dwn)<-react.tab0.grp.plot.title()
  gom.obj.dwn <- newGOM(genes.lists.dwn,genes.lists.dwn,
                        genome.size=gs.RNASeq)
  gom.obj.dwn
 })
 
 tab0.genes.lists <- reactive ({
  tab0.grp.name=react.tab0.grp.name()
  tab0.grp.plot.title=react.tab0.grp.plot.title()
  tab0.results=react.tab0.result()
  tab0.fdr.co <- react.tab0.fdr.cutoff()
  tab0.fc.cutoff <- react.tab0.fc.cutoff()
  tab0.grp.dwn <- tab0.grp.dwn ()
  tab0.grp.up <- tab0.grp.up ()
  # all regulated-genes
  tab0.grp=list()
  for(i in 1:tab0.ngroup){
   tab0.grp[[i]]=rbind(tab0.grp.up[[i]],tab0.grp.dwn[[i]])
  }
  
  tab0.genes.lists=lapply(tab0.grp,function(x) rownames(x))
  names(tab0.genes.lists)<-tab0.grp.plot.title
  tab0.genes.lists
 })
 
 tab0.gom.obj <- reactive({
  tab0.genes.lists <- tab0.genes.lists()
  tab0.gs.RNASeq <- tab0.gs.RNASeq()
  newGOM(tab0.genes.lists,tab0.genes.lists,genome.size=tab0.gs.RNASeq)
 })
 
 tab0.s4 <- reactive({
  # plot euler venn
  tab0.s4<-tab0.genes.lists()
  names(tab0.s4)<-react.tab0.grp.plot.title()
  tab0.s4
 })
 
 react.tab0.colors <- reactive({
  tab0.colors[1]=input$tab0.color.grp1
  tab0.colors[2]=input$tab0.color.grp2
  tab0.colors
 })
 tab0.s4.dwn <- reactive({
  # plot euler venn
  tab0.s4<-tab0.genes.lists.dwn()
  names(tab0.s4)<-react.tab0.grp.plot.title()
  tab0.s4
 })
 
 tab0.data.up <- reactive({
  gom.obj.up <- tab0.gom.obj.up()
  grp.plot.title <- react.tab0.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj.up,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[2,]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj.up,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[2,]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 tab0.data.dwn <- reactive({
  gom.obj.dwn <- tab0.gom.obj.dwn()
  grp.plot.title <- react.tab0.grp.plot.title()
  ## heatmap down regulation
  data <- getMatrix(gom.obj.dwn,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[2,]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj.dwn,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[2,]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 tab0.data <- reactive({
  gom.obj <- tab0.gom.obj()
  grp.plot.title <- react.tab0.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data[2,]
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp[2,]
  data$jaccard=tmp$value
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 output$tab0.plot1 <- renderPlot({
  venn.opts=react.tab0.venn.opts()
  colors=react.tab0.colors()
  cex=react.tab0.venn.cex()
  pad=react.tab0.venn.pad()
  set.seed(1)
  s4 <- tab0.s4.up()
  venn.up <- plot_euler(s4=s4,colors=colors,cex=cex,venn.opts=venn.opts,title="Up-regulation")
  s4 <- tab0.s4.dwn()
  venn.dwn <- plot_euler(s4=s4,colors=colors,cex=cex,venn.opts=venn.opts,title="Down-regulation")
  grid.arrange(venn.up,venn.dwn,padding=unit(pad,"line"),ncol=2,
               top="",left="",right="",bottom="")
 })
 
 output$tab0.plot2 <- renderPlot({
  data=tab0.data.up()
  hm.up <- heatmap_sigf_overlap(data)
  data=tab0.data.dwn()
  hm.dwn <- heatmap_sigf_overlap(data)
  pad=react.tab0.venn.pad()
  grid.arrange(hm.up,hm.dwn,padding=unit(pad,"line"),ncol=2,
               top="",left="",right="",bottom="")
 })
 
 output$tab0.plot3 <- renderPlot({	
  venn.opts=react.tab0.venn.opts()
  colors=react.tab0.colors()
  cex=react.tab0.venn.cex()
  pad=react.tab0.venn.pad()
  set.seed(1)
  s4 <- tab0.s4()
  plot_euler(s4=s4,colors=colors,cex=cex,venn.opts=venn.opts,title="All regulated genes")
 })
 
 output$tab0.plot4 <- renderPlot({
  data=tab0.data()
  heatmap_sigf_overlap(data)
 })	
 
 #### processing > 3 groups	
 
 react.tab3.grp.name <- reactive({
  grp.name <- sapply(grep("tab3\\.grp.+\\.name", x = names(input), value = TRUE),
                     function(x) input[[x]])
  grp.name=as.vector(grp.name[mixedorder(names(grp.name))])
  grp.name[grp.name!="NA"]
 })
 
 react.tab3.grp.plot.title <- reactive({
  grp.plot.title <- sapply(grep("tab3\\.grp.+\\.plot\\.title", x = names(input), value = TRUE),
                           function(x) input[[x]])
  grp.plot.title=as.vector(grp.plot.title[mixedorder(names(grp.plot.title))])
  grp.plot.title[grp.plot.title!=""]
 })
 
 react.tab3.color <- reactive({
  colors <- sapply(grep("tab3\\.color.+", x = names(input), value = TRUE),
                   function(x) input[[x]])
  colors=as.vector(colors[mixedorder(names(colors))])
  colors[colors!=""]
 })
 react.tab3.fdr.cutoff <- reactive({
  fdr.cutoff <- sapply(grep("tab3\\.fdr.+", x = names(input), value = TRUE),
                       function(x) input[[x]])
  fdr.cutoff=as.vector(fdr.cutoff[mixedorder(names(fdr.cutoff))])
  fdr.cutoff[fdr.cutoff!="NA"]
 })  
 react.tab3.fc.cutoff <- reactive({
  fc.cutoff <- sapply(grep("tab3\\.fc.+", x = names(input), value = TRUE),
                      function(x) input[[x]])
  fc.cutoff=as.vector(fc.cutoff[mixedorder(names(fc.cutoff))])
  fc.cutoff[fc.cutoff!="NA"]
 })
 react.tab3.meanDiff.cutoff <- reactive({
  meanDiff.cutoff <- sapply(grep("tab3\\.meanDiff.+", x = names(input), value = TRUE),
                            function(x) input[[x]])
  meanDiff.cutoff=as.vector(meanDiff.cutoff[mixedorder(names(meanDiff.cutoff))])
  meanDiff.cutoff[meanDiff.cutoff!="NA"]
 }) 
 react.tab3.venn.opts <- reactive({input$tab3.venn.opts})
 
 react.tab3.venn.cex <- reactive({input$tab3.venn.cex})
 
 react.tab3.venn.pad <- reactive({input$tab3.venn.pad})
 
 react.tab3.msigdb.species <- reactive({
  input$tab3.msigdb.species
 })
 
 react.tab3.nintersects.upset <- reactive({input$tab3.nintersects.upset})
 
 react.tab3.enrich.pval.co <- reactive({input$tab3.enrich.pval.co})
 
 react.tab3.chg.enrich.terms<- reactive({
  tmp=str_trim(unlist(strsplit(input$tab3.chg.enrich.terms,",")))
  tmp=strsplit(tmp,"=")
  shiny::req(!is.null(tmp))
  terms=sapply(tmp,"[[",2)
  names(terms)=sapply(tmp,"[[",1)
  dict = c("Mrna" = "mRNA", "Dna" = "DNA", "Rna" = "RNA", 
           "Trna" = "tRNA", "Mirna" = "miRNA", "Rrna" = "rRNA",
           "Atp" = "ATP","Adp" = "ADP","Snorna" = "snoRNA",
           "Lncrna" = "lncRNA","Snrna" = "snRNA",
           "Ncrna" = "ncRNA")
  
  if(length(terms)>0){
   dict = c(dict,terms)
  }
  dict
  })
 
 react.tab3.enrich.ncat <- reactive({
  input$tab3.enrich.ncat
 })
 
 react.tab3.enrich.plot.size <- reactive({
 input$tab3.enrich.plot.size
 })
 
 react.tab3.result <- reactive({
  grp.name=react.tab3.grp.name()
  projdir = react.setup.proj.dir()
  shiny::req(!grp.name=="")
  grp.plot.title=react.tab3.grp.plot.title()
  out.DESeq2 <- react.tab3.rdata()
  idx=match(grp.name,names(out.DESeq2$results))
  results=out.DESeq2$results[idx]
  # make sure grp list is generated before proceeding
  shiny::req(!is.null(rownames(results[[1]])[1]))
  # If rownames is ensebml id, change them to gene symbol
  if(sum(grepl("ENS",rownames(results[[1]])[1:10]))>2){
   in.file=paste0(gsub("_","",grp.name[1]),".complete.txt")
   df=read.table(paste0(projdir,"outputs/diff_analysis_rslt/tables/",in.file),
                sep="\t",header=TRUE)
   gene.symbol=df$geneSymbol
   # for genes with same gene symbol but different id, keep the first duplicate
   # change the rest to id
   dup <- duplicated(gene.symbol, fromLast = TRUE)
   gene.symbol[dup] <- NA
   # for genes without gene symbol, set id as gene symbol
   gene.symbol[is.na(gene.symbol)]=df$Id[is.na(gene.symbol)]
   results=lapply(results,function(x) {
    rownames(x)=gene.symbol 
    return(x)
    })
  }
  results
 })
 
 react.tab3.normCts <- reactive({
  out.DESeq2 <- react.tab3.rdata()
  normCts=counts(out.DESeq2$dds,normalized=TRUE)
 })
 
 tab3.gs.RNASeq <- reactive({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  gs.RNASeq <- max(sapply(results,function(x) sum(!is.na(x$padj))))
 })
 
 tab3.grp.up <- reactive ({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  normCts=react.tab3.normCts()
  fdr.co <- react.tab3.fdr.cutoff()
  fc.cutoff <- react.tab3.fc.cutoff()
  meanDiff.cutoff <- react.tab3.meanDiff.cutoff()
  out.DESeq2=react.tab3.rdata()
  
  # needed to correctly differentiate sample names
  # i.e. iEF vs iEFempty
  rep.name=levels(colData(out.DESeq2$dds)$rep)
  # if rep is not used in DESeq, it will remain a character
  # therefore rep.name will be NULL
  # so instead get the unique rep
  if(is.null(rep.name)) rep.name=unique(colData(out.DESeq2$dds)$rep)
  
  # up-regulated genes
  grp.up=list()
  for(i in 1:length(grp.name)){
   
   # filter out by mean difference
   grps=strsplit(grp.name[i],"_vs_")
   treat.grp=grps[[1]][1]
   ref.grp=grps[[1]][2]
   # fixed bug: below code will get both UASVN2_1055A and 1055A samples while
   # grepping only 1055A 
   #treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
   #ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
   # adding rep.name so it will get the correct samples.
   treat.mean=rowMeans(normCts[,which(colnames(normCts) %in% paste0(treat.grp,rep.name))])
   ref.mean=rowMeans(normCts[,which(colnames(normCts) %in% paste0(ref.grp,rep.name))])
   diff.mean=treat.mean-ref.mean
   # grp.up[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
   # results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) &
   # !is.na(results[[grp.name[i]]]$padj)& diff.mean >= meanDiff.cutoff[i],]
   row.pass=which(results[[grp.name[i]]]$padj<=fdr.co[i] &
                   results[[grp.name[i]]]$log2FoldChange>(log2(fc.cutoff[i])) &
                   !is.na(results[[grp.name[i]]]$padj)& diff.mean >= meanDiff.cutoff[i])
   grp.up[[i]]=results[[grp.name[i]]][row.pass,]
   grp.up[[i]]$row.pass = NULL
   grp.up[[i]]$norm.diff.mean = NULL
   if(length(row.pass)>0){
    grp.up[[i]]$row.pass = row.pass
    grp.up[[i]]$norm.diff.mean = diff.mean[row.pass]
   }
  }
  grp.up
 })
 
 tab3.grp.dwn <- reactive({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  normCts = react.tab3.normCts()
  fdr.co <- react.tab3.fdr.cutoff()
  fc.cutoff <- react.tab3.fc.cutoff()
  meanDiff.cutoff <- react.tab3.meanDiff.cutoff()
  out.DESeq2=react.tab3.rdata()
  
  # get the replicate names
  # needed to correctly differentiate sample names
  # i.e. iEF vs iEFempty
  rep.name=levels(colData(out.DESeq2$dds)$rep)
  # if rep is not used in DESeq, it will remain a character
  # therefore rep.name will be NULL
  # so instead get the unique rep
  if(is.null(rep.name)) rep.name=unique(colData(out.DESeq2$dds)$rep)
  
  # down-regulated genes
  grp.dwn=list()
  for(i in 1:length(grp.name)){
   # filter out by mean difference
   grps=strsplit(grp.name[i],"_vs_")
   treat.grp=grps[[1]][1]
   ref.grp=grps[[1]][2]
   # fixed bug: below code will get both UASVN2_1055A and 1055A samples while
   # grepping only 1055A 
   #treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
   #ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
   # adding rep.name so it will get the correct samples.
   treat.mean=rowMeans(normCts[,which(colnames(normCts) %in% paste0(treat.grp,rep.name))])
   ref.mean=rowMeans(normCts[,which(colnames(normCts) %in% paste0(ref.grp,rep.name))])
   diff.mean=ref.mean-treat.mean
   grp.dwn[[i]]=results[[grp.name[i]]][results[[grp.name[i]]]$padj<=fdr.co[i] &
                                        results[[grp.name[i]]]$log2FoldChange<(-log2(fc.cutoff[i])) &
                                        !is.na(results[[grp.name[i]]]$padj) & diff.mean >= meanDiff.cutoff[i],]
   row.pass=which(results[[grp.name[i]]]$padj<=fdr.co[i] &
                   results[[grp.name[i]]]$log2FoldChange<(-log2(fc.cutoff[i])) &
                   !is.na(results[[grp.name[i]]]$padj) & diff.mean >= meanDiff.cutoff[i])
   grp.dwn[[i]]=results[[grp.name[i]]][row.pass,]
   grp.dwn[[i]]$row.pass = NULL
   grp.dwn[[i]]$norm.diff.mean = NULL
   if(length(row.pass)>0){
    grp.dwn[[i]]$row.pass = row.pass
    grp.dwn[[i]]$norm.diff.mean = diff.mean[row.pass]
   }
  }#for(i in 1:length(grp.name))
  print("length(grp.dwn)")
  print(length(grp.dwn))
  grp.dwn
 })
 
 tab3.genes.lists <- reactive ({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  grp.dwn <- tab3.grp.dwn ()
  grp.up <- tab3.grp.up ()
  
  # all regulated-genes
  grp=list()
  for(i in 1:length(grp.name)){
   grp[[i]]=rbind(grp.up[[i]],grp.dwn[[i]])
  }
  
  genes.lists=lapply(grp,function(x) rownames(x))
  names(genes.lists)<-grp.plot.title
  genes.lists
 })
 
 tab3.genes.lists.up <- reactive ({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  grp.up <- tab3.grp.up ()
  genes.lists.up=lapply(grp.up,function(x) rownames(x))
  names(genes.lists.up) <- grp.plot.title
  genes.lists.up
 })
 
 tab3.genes.lists.dwn <- reactive ({
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  grp.dwn <- tab3.grp.dwn ()
  genes.lists.dwn=lapply(grp.dwn,function(x) rownames(x))
  names(genes.lists.dwn) <- grp.plot.title
  genes.lists.dwn
 })
 
 tab3.gom.obj <- reactive({
  genes.lists <- tab3.genes.lists()
  gs.RNASeq <- tab3.gs.RNASeq()
  
  newGOM(genes.lists,genes.lists,genome.size=gs.RNASeq)
 })
 
 tab3.gom.obj.up <- reactive({
  genes.lists <- tab3.genes.lists.up()
  gs.RNASeq <- tab3.gs.RNASeq()
  newGOM(genes.lists,genes.lists,genome.size=gs.RNASeq)
 })
 
 tab3.gom.obj.dwn <- reactive({
  genes.lists <- tab3.genes.lists.dwn()
  gs.RNASeq <- tab3.gs.RNASeq()
  newGOM(genes.lists,genes.lists,genome.size=gs.RNASeq)
 })
 
 tab3.data <- reactive({
  gom.obj <- tab3.gom.obj()
  grp.plot.title <- react.tab3.grp.plot.title()
  
  # calculate overlap
  data <- getMatrix(gom.obj,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp
  data$jaccard=tmp$value
  # categorizing pvalues
  data$pval.cat=ifelse(data$pval==0,pval.names[1],
               ifelse(data$pval<=0.05,pval.names[2],
                       ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 tab3.data.up <- reactive({
  gom.obj <- tab3.gom.obj.up()
  grp.plot.title <- react.tab3.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp
  data$jaccard=tmp$value
  # categorizing pvalues
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 tab3.data.dwn <- reactive({
  gom.obj <- tab3.gom.obj.dwn()
  grp.plot.title <- react.tab3.grp.plot.title()
  # calculate overlap
  data <- getMatrix(gom.obj,"pval")
  data <- melt(data,id=grp.plot.title)
  data <- data
  names(data)[names(data)=="value"]="pval"
  tmp <- getMatrix(gom.obj,"Jaccard")
  tmp <- melt(tmp,id=grp.plot.title)
  tmp <- tmp
  data$jaccard=tmp$value
  # categorizing pvalues
  data$pval.cat=ifelse(data$pval==0,pval.names[1],ifelse(data$pval<=0.05,pval.names[2],
                                                         ifelse(data$pval<0.5,pval.names[3],pval.names[4])))
  data$pval.cat=factor(data$pval.cat,levels=pval.names)
  data
 })
 
 tab3.s4 <- reactive({
  # input for euler venn
  s4 <-tab3.genes.lists()
  grp.plot.title <-react.tab3.grp.plot.title()
  
  names(s4) <- grp.plot.title
  s4
 })
 
 tab3.s4.up <- reactive({
  # input for euler venn
  s4 <-tab3.genes.lists.up()
  grp.plot.title <-react.tab3.grp.plot.title()
  
  names(s4) <- grp.plot.title
  s4
 })
 
 tab3.s4.dwn <- reactive({
  # input for euler venn
  s4 <-tab3.genes.lists.dwn()
  grp.plot.title <-react.tab3.grp.plot.title()
  names(s4) <- grp.plot.title
  s4
 })
 
 tab3.compare.df <- reactive({
  # input for clusterProfiler
  grp = tab3.grp.up()
  grp.name = react.tab3.grp.name()
  grp.plot.title =react.tab3.grp.plot.title()
  compare.df = NULL
  for(i in 1:length(grp.plot.title)){
   compare.df=rbind(compare.df,data.frame(SYMBOL=rownames(grp[[i]]),
                                          group1=rep(grp.plot.title[i],dim(grp[[i]])[1]),
                                          group2=rep("Up-regulated",dim(grp[[i]])[1])))
  }
  grp = tab3.grp.dwn()
  for(i in 1:length(grp.plot.title)){
   compare.df=rbind(compare.df,data.frame(SYMBOL=rownames(grp[[i]]),
                                          group1=rep(grp.plot.title[i],dim(grp[[i]])[1]),
                                          group2=rep("Down-regulated",dim(grp[[i]])[1])))
  }
  compare.df
 })
 
  # preprocess multiple search term
 react.search.term <- reactive({
  if(!is.null(input$tab3.search.term) && !input$tab3.search.term =="")
   paste(str_trim(unlist(strsplit(input$tab3.search.term,","))),collapse="|")
 })
 
 # Filter data based on the search input
 react.tab3.filtered.data <- reactive({
  search_term <- react.search.term()
  data_list <- react.tab3.expr.tbl()
  current_page <- input$currentPage
  if (is.null(search_term) || search_term == "") {
   rows_per_page <- 10
   start_row <- (current_page-1)*rows_per_page+1
   end_row <- min(current_page*rows_per_page,sum(sapply(data_list, nrow)))
   sliced_data_lst <- lapply(data_list, function(df) {
     slice(df, start_row:end_row)
    })
   return(sliced_data_lst)
  } else {
   filtered_list <- lapply(data_list, function(df) {
    df[grep(search_term, df$Genes, ignore.case = TRUE), , drop = FALSE]
   })
   return(filtered_list)
  }
 })

 # Dynamically generate tables with titles
 output$tables <- renderUI({
  data_list <- react.tab3.expr.tbl()
  lapply(names(data_list), function(tab_name) {
   table_id <- paste0("table_", tab_name)
   tagList(
    h4(paste("Gene expression for", tab_name)),
    tableOutput(table_id),
    br()
   )
  })
 })
 
 # Render the filtered data tables
 observe({
  data_list <- react.tab3.expr.tbl()
  filtered.data <- react.tab3.filtered.data()
  for (tab_name in names(data_list)) {
   # Need local so that each item gets its own number. Without it, the value
   # of tab_name in the renderTable() will be the same across all instances, because
   # of when the expression is evaluated.
   # https://gist.github.com/wch/5436415/
   local({
    tab_name <- tab_name
   table_id <- paste0("table_", tab_name)
   output[[table_id]] <- renderTable({
     as.data.frame(filtered.data[[tab_name]])
   })
   })#local
  }
 })

 # Insert the right number of plot output objects into the web page
 output$tab3.volcano.plots <- renderUI({
  grp.name=react.tab3.grp.name()
  genes.lists.dwn = tab3.genes.lists.dwn()
  print("genes.lists.dwn")
  print(length(genes.lists.dwn))
  genes.lists.up = tab3.genes.lists.up()
  plot_output_list<-list()
  for (i in 1:length(grp.name)){
   plot_output_list[[i]] <- plotOutput(paste0("tab3.volcano.grp",i),height="600px")
  }
  # Convert the list to a tagList - this is necessary for the list of items
  # to display properly.
  do.call(tagList, plot_output_list)
 })
 
 # Insert the right number of plot output objects into the web page
 output$tab3.venn.plots <- renderUI({
  grp.name=react.tab3.grp.name()
  genes.lists.dwn = tab3.genes.lists.dwn()
  genes.lists.up = tab3.genes.lists.up()
  plot_output_list<-list()
  plot.i=1
  if(length(grp.name)<8){
   if(sum(unlist(sapply(genes.lists.up,function(x)length(x))))>0){
    plot_output_list[[length(grp.name)+plot.i]] <- plotOutput("tab3.venn.up")
    plot.i = plot.i +1
   }
   if(sum(unlist(sapply(genes.lists.dwn,function(x)length(x))))>0){
    plot_output_list[[length(grp.name)+plot.i]] <- plotOutput("tab3.venn.dwn")
    plot.i=plot.i+1
   }
   if(length(grp.name)>1)
    plot_output_list[[length(grp.name)+plot.i]] <- plotOutput("tab3.hm")
  }else{
   plot_output_list[[length(grp.name)+plot.i]] <- plotOutput("tab3.hm")
  }
  
  # Convert the list to a tagList - this is necessary for the list of items
  # to display properly.
  do.call(tagList, plot_output_list)
 })
 
 # # Call renderPlot for each one. Plots are only actually generated when they
 # # are visible on the web page.
 
 react.tab3.hilite.genes <- reactive({
  unlist(sapply(strsplit(input$tab3.hilite.genes,","),function(x)str_trim(x)))
 })
 
 react.tab3.expr.tbl <- reactive({
  shiny::req(input$tabset=="Plots")
  grp.name=react.tab3.grp.name()
  grp.plot.title=react.tab3.grp.plot.title()
  results=react.tab3.result()
  normCts=react.tab3.normCts()
  fdr.co <- react.tab3.fdr.cutoff()
  fc.cutoff <- react.tab3.fc.cutoff()
  meanDiff.cutoff <- react.tab3.meanDiff.cutoff()
  out.DESeq2=react.tab3.rdata()
  
  # get the replicate names
  # needed to correctly differentiate sample names
  # i.e. iEF vs iEFempty
  rep.name=levels(colData(out.DESeq2$dds)$rep)
  # if rep is not used in DESeq, it will remain a character
  # therefore rep.name will be NULL
  # so instead get the unique rep
  if(is.null(rep.name)) rep.name=unique(colData(out.DESeq2$dds)$rep)
  
  expr.tbl <- list()
  for(i in 1:length(grp.name)){
   # get indiv grpname
   grps=strsplit(grp.name[i],"_vs_")
   treat.grp=grps[[1]][1]
   ref.grp=grps[[1]][2]
   print("treat.grp")
   print(treat.grp)
   # fixed bug: below code will get both UASVN2_1055A and 1055A samples while
   # grepping only 1055A 
   #treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
   #ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
   # adding rep.name so it will get the correct samples.
   treat.mean=rowMeans(normCts[,which(colnames(normCts) %in% paste0(treat.grp,rep.name))])
   ref.mean=rowMeans(normCts[,which(colnames(normCts) %in% paste0(ref.grp,rep.name))])
   diff.mean=treat.mean-ref.mean
   print("colnames(normCts)")
   print(colnames(normCts))
   print("paste0(treat.grp,rep.name)")
   print(paste0(treat.grp,rep.name))
   print("paste0(ref.grp,rep.name)")
   print(paste0(ref.grp,rep.name))
   print("levels(colData(out.DESeq2$dds)$rep)")
   print(levels(colData(out.DESeq2$dds)$rep))
   data=data.frame(
    rownames(results[[grp.name[i]]]),
    results[[grp.name[i]]]$log2FoldChange,
    results[[grp.name[i]]]$padj,
    diff.mean,
    ref.mean,
    treat.mean)
   names(data)=c("Genes","log2FC","FDR","difference",ref.grp,
                 treat.grp)
   data <- data %>% 
    mutate(
     Expression = 
      case_when(log2FC >= log2(fc.cutoff[i]) & 
                 FDR <= fdr.co[i] & 
                 diff.mean >= meanDiff.cutoff[i] ~ "Up-regulated",
                log2FC < -log2(fc.cutoff[i]) &
                 FDR <= fdr.co[i] & 
                 diff.mean <= -meanDiff.cutoff[i] ~ "Down-regulated",
                TRUE ~ "NS")
    ) %>%
    arrange(round(FDR,digits=3),desc(abs(log2FC)))
   expr.tbl[[i]] <- data
  }
  names(expr.tbl)=grp.plot.title
  expr.tbl
 })
 
 observe({
 grp.name=react.tab3.grp.name()
 grp.plot.title=react.tab3.grp.plot.title()
 results=react.tab3.result()
 normCts=react.tab3.normCts()
 fdr.co <- react.tab3.fdr.cutoff()
 fc.cutoff <- react.tab3.fc.cutoff()
 meanDiff.cutoff <- react.tab3.meanDiff.cutoff()
 hilite.genes <- react.tab3.hilite.genes()
 out.DESeq2=react.tab3.rdata()
 
 # get the replicate names
 # needed to correctly differentiate sample names
 # i.e. iEF vs iEFempty
 rep.name=levels(colData(out.DESeq2$dds)$rep)
 # if rep is not used in DESeq, it will remain a character
 # therefore rep.name will be NULL
 # so instead get the unique rep
 if(is.null(rep.name)) rep.name=unique(colData(out.DESeq2$dds)$rep)
 
 for(i in 1:length(grp.name)){
  # Need local so that each item gets its own number. Without it, the value
  # of i in the renderPlot() will be the same across all instances, because
  # of when the expression is evaluated.
  # https://gist.github.com/wch/5436415/
  local({
   my_i <- i
 output[[paste0("tab3.volcano.grp",my_i)]] <- renderPlot({
   # get indiv group name
   grps=strsplit(grp.name[my_i],"_vs_")
   treat.grp=grps[[1]][1]
   ref.grp=grps[[1]][2]
   # fixed bug: below code will get both UASVN2_1055A and 1055A samples while
   # grepping only 1055A 
   #treat.mean=rowMeans(normCts[,grep(treat.grp,colnames(normCts))])
   #ref.mean=rowMeans(normCts[,grep(ref.grp,colnames(normCts))])
   # adding rep.name so it will get the correct samples.
   treat.mean=rowMeans(normCts[,which(colnames(normCts) %in% paste0(treat.grp,rep.name))])
   ref.mean=rowMeans(normCts[,which(colnames(normCts) %in% paste0(ref.grp,rep.name))])
   diff.mean=treat.mean-ref.mean
   data=data.frame(
    Genes=rownames(results[[grp.name[my_i]]]),
    logFC=results[[grp.name[my_i]]]$log2FoldChange,
    FDR=results[[grp.name[my_i]]]$padj,
    diff.mean=diff.mean)
   #print("data$Genes")
   #print(data[data$Genes =="polr3gla",])
   data <- data %>% 
    mutate(
     Expression = 
      case_when(logFC >= log2(fc.cutoff[my_i]) & 
                 FDR <= fdr.co[my_i] & 
                 diff.mean >= meanDiff.cutoff[my_i] ~ "Up-regulated",
                logFC < -log2(fc.cutoff[my_i]) &
                 FDR <= fdr.co[my_i] & 
                 diff.mean <= -meanDiff.cutoff[my_i] ~ "Down-regulated",
                TRUE ~ "NS")
    )
   #print("Down-regulated")
   #print(data[data$Expression =="Down-regulated",])
   # making sure colors are correctly assigned to expression
   data$Expression <- factor(data$Expression,levels=c("NS","Up-regulated","Down-regulated"))
   total_up=sum(data$Expression == "Up-regulated")
   total_down=sum(data$Expression == "Down-regulated")
   # change FDR=0 so it can be graphed correctly
   # remove FDR=na
   data=data[!is.na(data$FDR),]
   data$logFDR=-log(data$FDR+min(c(data$FDR[data$FDR>0],1e-32)),10)
   #print("Down-regulated")
   #print(data[data$Expression =="Down-regulated",])
   genes_show = hilite.genes
   genes_show_data <- dplyr::bind_rows(
    data %>%
     filter(Genes %in% genes_show)
   )
   xlim=c(-max(abs(data$logFC))*1.1,max(abs(data$logFC))*1.1)
   ylim=c(0,max(data$logFDR,na.rm=TRUE)*1.1)
   p <- ggplot(data, aes(logFC, logFDR)) +
    geom_hline(yintercept=-log10(fdr.co[my_i]), col="#5d5d5d",linetype="dashed")+
    geom_point(aes(color = Expression), size = 4,alpha=4/5) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"FDR")) +
    labs(
     subtitle=paste("Up-regulated:",total_up, " | Down-regulated:", total_down)
    )+
    scale_color_manual(values = c("NS"="gray50","Down-regulated"=tab3.colors[2],"Up-regulated"=tab3.colors[6])) +
    guides(colour = guide_legend(override.aes = list(size=2))) + 
    geom_label_repel(data = genes_show_data,mapping=aes(logFC,logFDR,label=Genes),
                     size=5,min.segment.length = 0) +
    theme_classic(base_size = 18) +
    theme(panel.grid.major=element_line(color = "#EBEBEB"),
          panel.grid.minor=element_line(color = "#EBEBEB"),
          # margin add space below title
          plot.title=element_text(hjust=0.5,margin=margin(b=20)))+
    ggtitle(paste0("Volcano plot for ",grp.plot.title[my_i])) +
    scale_x_continuous(limits=xlim)
   if(fc.cutoff[my_i]>1){
    p <- p + 
     geom_vline(xintercept=c(-log2(fc.cutoff[my_i]), log2(fc.cutoff[my_i])), col="#5d5d5d",linetype="dashed")
   }
   vals.plot.volcano[[paste0("tab3.volcano.grp",my_i)]] = p
   vals.plot.volcano[[paste0("tab3.volcano.grp",my_i)]]
 })#renderPlot  
  })#local
  }#for(i in 1:length(grp.name))
 })#observe

 output[["tab3.venn.up"]] <- renderPlot({
  venn.opts=react.tab3.venn.opts()
  grp.name=react.tab3.grp.name()
  colors=react.tab3.color()
  cex=react.tab3.venn.cex()
  pad=react.tab3.venn.pad()
  set.seed(1)
  s4 <- tab3.s4.up()
  # if all lists contains 0 genes then do not proceeds
  nonzero_s4=sum(sapply(s4,function(x) length(x)!=0))
  shiny::req(nonzero_s4>0)
  
  # remove label if requested
  # set names of s4 for venn library
  # note: snames is ignored when s4 has names
  save_s4_names=names(s4)
    if('Labels' %in% venn.opts){
     names(s4)=str_wrap(gsub("_"," ",names(s4)),10)
     ilabels="counts"
     labels=list(names(s4),cex=cex)
  }
  else{
   labels=""
   names(s4)=sapply(1:length(s4),function(x) strrep(" ",x))
   ilabels="counts"
  }
  if(!'Numbers' %in% venn.opts) ilabels=NULL
  
  vals.plot$venn.up1<-venn(s4,zcolor=colors,opacity=.8,box=FALSE,
                           ilcs = 0.8, sncs = 1,ggplot=TRUE,ilabels=ilabels)
  # only for euler
  if(!'Labels' %in% venn.opts && 'Legend' %in% venn.opts){
   names(s4)=save_s4_names
   labels=""
  }
  vals.plot$venn.up2<-plot_euler(s4=s4,colors=colors,cex=cex,
                                 venn.opts=venn.opts,title="",labels=labels)

  grid.arrange(vals.plot$venn.up1,vals.plot$venn.up2,ncol=2,
               padding=unit(pad,"line"),
               top=textGrob("Up-regulated Genes",gp=gpar(fontsize=20)))
 })
 
 output[["tab3.venn.dwn"]] <- renderPlot({
  venn.opts=react.tab3.venn.opts()
  grp.name=react.tab3.grp.name()
  colors=react.tab3.color()
  cex=react.tab3.venn.cex()
  pad=react.tab3.venn.pad()
  set.seed(1)
  s4 <- tab3.s4.dwn()
  # if all lists contains 0 genes then do not proceeds
  nonzero_s4=sum(sapply(s4,function(x) length(x)!=0))
  shiny::req(nonzero_s4>0)
  # remove label if requested
  # set names of s4 for venn library
  # note: snames is ignored when s4 has names
  save_s4_names=names(s4)
  if('Labels' %in% venn.opts){
   names(s4)=str_wrap(gsub("_"," ",names(s4)),10)
   ilabels="counts"
   labels=list(names(s4),cex=cex)
  }
  else{
   labels=""
   names(s4)=sapply(1:length(s4),function(x) strrep(" ",x))
   ilabels="counts"
  }
  if(!'Numbers' %in% venn.opts) ilabels=NULL
  vals.plot$venn.dwn1 <- venn(s4,zcolor=colors,opacity=.8,box=FALSE,ilcs = 0.8,
                              sncs = 1,ggplot=TRUE,ilabels=ilabels)
  # only for euler
  if(!'Labels' %in% venn.opts && 'Legend' %in% venn.opts){
   names(s4)=save_s4_names
   labels=""
  }
  vals.plot$venn.dwn2<-plot_euler(s4=s4,colors=colors,cex=cex,
                                  venn.opts=venn.opts,title="",labels=labels)



  grid.arrange(vals.plot$venn.dwn1,vals.plot$venn.dwn2,ncol=2,
               padding=unit(pad,"line"),
               top=textGrob("Down-regulated Genes",gp=gpar(fontsize=20)))
 })
 
 output[["tab3.hm"]] <- renderPlot({
  pad = react.venn.pad()
  data = tab3.data.up()
  grp.name = react.tab3.grp.name()
  shiny::req(length(grp.name)>1)
  vals.plot$hm.up <- heatmap_sigf_overlap(data,"Up-regulated genes")
  data=tab3.data.dwn()
  vals.plot$hm.dwn<- heatmap_sigf_overlap(data,"Down-regulated genes")
  grid.arrange(vals.plot$hm.dwn,vals.plot$hm.up,ncol=2)
 })
 
 output$tab3.plot1 <- renderPlot({
  s4 <- tab3.s4.up()
  grp.name <- react.tab3.grp.name()
  # only draw upset up-regulated if there are more than 1 group with non-empty list of up genes
  shiny::req(length(grp.name)>1 && sum(sapply(s4,function(x)length(x)!=0))>1)
  nintersects <- react.tab3.nintersects.upset()
  vals.plot$upset.up <- as.ggplot(upset(fromList(s4), order.by = "freq",
                                        nsets=length(grp.name), nintersects=nintersects, 
                                        mainbar.y.label = "# Up-regulated Genes", sets.x.label = "Significant Genes", 
                                        point.size=2.8,line.size=1,
                                        #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
                                        text.scale = c(1.8, 1.8, 1.4, 1.8, 1.8, 1.5)))
  vals.plot$upset.up
 })

 
 output$tab3.plot2 <- renderPlot({
  s4 <- tab3.s4.dwn()
  grp.name <- react.tab3.grp.name()
  # only draw upset if there are more than two groups with non-empty list of diff genes
  shiny::req(length(grp.name)>1 && sum(sapply(s4,function(x)length(x)!=0))>1)
  nintersects <- react.tab3.nintersects.upset()
  vals.plot$upset.dwn <- as.ggplot(upset(fromList(s4), order.by = "freq",nsets=length(grp.name), 
                                         nintersects=nintersects,
                                         mainbar.y.label = "# Down-regulated Genes", sets.x.label = "Significant Genes", 
                                         point.size=2.8,line.size=1,
                                         #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
                                         text.scale = c(1.8, 1.8, 1.4, 1.8, 1.8, 1.5)))
  vals.plot$upset.dwn
 })
 
 # gene sets from msigdb
 # # oncogenic signature gene
 # m_t2g.c6 <- msigdbr(species = msigdb.species, category = "C6") %>%
 #  dplyr::select(gs_name, gene_symbol)
 # # hallmark gene sets
 # m_t2g.h <- msigdbr(species = msigdb.species, category = "H") %>%
 #  dplyr::select(gs_name, gene_symbol)
 # m_t2g.c2.kegg <- msigdbr(species = msigdb.species, category = "C2",subcategory = "CP:KEGG") %>%
 #  dplyr::select(gs_name, gene_symbol)
 # m_t2g.c2.reactome <- msigdbr(species = msigdb.species, category = "C2",subcategory = "CP:REACTOME") %>%
 #  dplyr::select(gs_name, gene_symbol)
 # m_t2g.c2.biocarta <- msigdbr(species = msigdb.species, category = "C2",subcategory = "CP:BIOCARTA") %>%
 #  dplyr::select(gs_name, gene_symbol)
 
 
  output$tab3.plot3.1 <- renderPlot({
   shiny::req(input$gen.go)
   isolate({
   withProgress(message="Generating enrichment plots",{
    enrich.pval.co <- react.tab3.enrich.pval.co()
    compare.df <- tab3.compare.df()  
    grp.plot.title <- react.tab3.grp.plot.title()
    chg.enrich.terms <- react.tab3.chg.enrich.terms()
    enrich.ncat <- react.tab3.enrich.ncat()
    # Using clusterProfiler to perform hypergeometric test on msigdb signatures
    msigdb.species <- react.tab3.msigdb.species()
    msig.gene.set <- msigdbr(species = msigdb.species, category = "C5",subcategory = "MF") %>%
     dplyr::select(gs_name, gene_symbol)
    msig.name ="MSigDB GO Molecular Function"

    formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
                                  TERM2GENE=msig.gene.set,pvalueCutoff=enrich.pval.co,
                                  pAdjustMethod="BH")

    # re-arrange datasets using factor
    # and do pathway analysis using up- and down-regulated genes separately
    vals.plot$msig.mf <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name,
                                 showCategory = enrich.ncat) + 
     facet_grid(~group2) +
     scale_y_discrete(labels=function(x) 
      str_wrap(str_replace_all(str_to_title(tolower(gsub("_"," ",gsub("GOMF_","",x)))),
                               chg.enrich.terms), width=40)) +
     scale_x_discrete(labels=function(x) 
      str_wrap(gsub("_"," ",x),width=10)) +
     #scale_color_distiller(palette = 'Blues') +
     scale_color_gradientn(colors=rev(pal(10)),
                           limits=c(0,enrich.pval.co))+
     theme(strip.text=element_text(size=14))
    vals.plot$msig.mf
   })#withProgress
   }) #isolate
  },height=function(){100*react.tab3.enrich.plot.size()}) #renderPlot
 
  output$tab3.plot3.2 <- renderPlot({
   shiny::req(input$gen.go)
   shiny::isolate({
   withProgress(message="Generating enrichment plots",{
    enrich.pval.co <- isolate(react.tab3.enrich.pval.co())
    compare.df <- isolate(tab3.compare.df())
    grp.plot.title <- isolate(react.tab3.grp.plot.title())
    # Using clusterProfiler to perform hypergeometric test on msigdb signatures
    msigdb.species <- isolate(react.tab3.msigdb.species())
    msig.gene.set <- msigdbr(species = msigdb.species, category = "C5",subcategory = "BP") %>%
     dplyr::select(gs_name, gene_symbol)
    msig.name ="MSigDB GO Biological Process"
    chg.enrich.terms <- react.tab3.chg.enrich.terms()
    enrich.ncat <- react.tab3.enrich.ncat()
    formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
                                  TERM2GENE=msig.gene.set,pvalueCutoff=enrich.pval.co,
                                  pAdjustMethod="BH")
    vals.plot$msig.bp <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name,
                                 showCategory = enrich.ncat) +
     facet_grid(~group2) +
     scale_y_discrete(labels=function(x)
      str_wrap(str_replace_all(str_to_title(tolower(gsub("_"," ",gsub("GOBP_","",x)))),chg.enrich.terms), width=40)) +
     scale_x_discrete(labels=function(x) str_wrap(gsub("_"," ",x),width=10)) +
     scale_color_gradientn(colors=rev(pal(10)),
                           limits=c(0,enrich.pval.co))+
     theme(strip.text=element_text(size=14))
    vals.plot$msig.bp
   })#withProgress
   }) #isolate
  },height=function(){100*react.tab3.enrich.plot.size()})

  output$tab3.plot3.3 <- renderPlot({
   shiny::req(input$gen.go)
   shiny::isolate({
   withProgress(message="Generating enrichment plots",{
    enrich.pval.co <- isolate(react.tab3.enrich.pval.co())
    compare.df <- isolate(tab3.compare.df())
    grp.plot.title <- isolate(react.tab3.grp.plot.title())
    chg.enrich.terms <- isolate(react.tab3.chg.enrich.terms())
    enrich.ncat <- react.tab3.enrich.ncat()
    # Using clusterProfiler to perform hypergeometric test on msigdb signatures
    msigdb.species <- isolate(react.tab3.msigdb.species())
    msig.gene.set <- msigdbr(species = msigdb.species, category = "C5",subcategory = "CC") %>%
     dplyr::select(gs_name, gene_symbol)
    msig.name ="MSigDB GO Cellular Component"
    formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
                                  TERM2GENE=msig.gene.set,pvalueCutoff=enrich.pval.co,
                                  pAdjustMethod="BH")
    vals.plot$msig.cc <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name,
                                showCategory=enrich.ncat) +
     facet_grid(~group2) +
     scale_y_discrete(labels=function(x)
      str_wrap(str_replace_all(str_to_title(tolower(gsub("_"," ",gsub("GOCC_","",x)))),chg.enrich.terms), width=40)) +
     scale_x_discrete(labels=function(x) str_wrap(gsub("_"," ",x),width=10)) +
     scale_color_gradientn(colors=rev(pal(10)),
                           limits=c(0,enrich.pval.co))+
     theme(strip.text=element_text(size=14))
    vals.plot$msig.cc
   }) #withProgress
   }) #isolate
  },height=function(){100*react.tab3.enrich.plot.size()})


  output$tab3.plot3.4 <- renderPlot({
   shiny::req(input$gen.go)
   shiny::isolate({
    withProgress(message="Generating enrichment plots",{
    enrich.pval.co <- isolate(react.tab3.enrich.pval.co())
    compare.df <- isolate(tab3.compare.df())
    grp.plot.title <- isolate(react.tab3.grp.plot.title())
    chg.enrich.terms <- isolate(react.tab3.chg.enrich.terms())
    enrich.ncat <- react.tab3.enrich.ncat()
    # Using clusterProfiler to perform hypergeometric test on msigdb signatures
    msigdb.species <- isolate(react.tab3.msigdb.species())
    msig.gene.set <- msigdbr(species = msigdb.species, category = "C2") %>%
     dplyr::select(gs_name, gene_symbol)
    msig.name ="MSigDB Curated Gene Sets"
    formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
                                  TERM2GENE=msig.gene.set,
                                  pvalueCutoff=enrich.pval.co,pAdjustMethod="BH")
    vals.plot$msig.curate <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name,
                                     showCategory=enrich.ncat) + facet_grid(~group2) +
     scale_y_discrete(labels=function(x)
            str_wrap(str_replace_all(str_to_title(tolower(gsub("_"," ",x))),chg.enrich.terms), width=40)) +
     scale_x_discrete(labels=function(x) str_wrap(gsub("_"," ",x),width=10)) +
     scale_color_gradientn(colors=rev(pal(10)),
                           limits=c(0,enrich.pval.co))+
     theme(strip.text=element_text(size=14))
    vals.plot$msig.curate
   }) # withProgress
   }) #isolate
  },height=function(){100*react.tab3.enrich.plot.size()})
 
  output$tab3.plot3.5 <- renderPlot({
   shiny::req(input$gen.go)
   shiny::isolate({
    withProgress(message="Generating enrichment plots",{
     enrich.pval.co <- isolate(react.tab3.enrich.pval.co())
     compare.df <- isolate(tab3.compare.df())
     grp.plot.title <- isolate(react.tab3.grp.plot.title())
     chg.enrich.terms <- isolate(react.tab3.chg.enrich.terms())
     enrich.ncat <- react.tab3.enrich.ncat()
     # Using clusterProfiler to perform hypergeometric test on msigdb signatures
     msigdb.species <- isolate(react.tab3.msigdb.species())
     msig.gene.set <- msigdbr(species = msigdb.species, category = "H") %>%
      dplyr::select(gs_name, gene_symbol)
     msig.name ="MSigDB Hallmark Gene Sets"
     formula_res <- compareCluster(SYMBOL~group1+group2, data=compare.df, fun="enricher",
                                   TERM2GENE=msig.gene.set,
                                   pvalueCutoff=enrich.pval.co,pAdjustMethod="BH")
     vals.plot$msig.h <- dotplot(formula_res,x=~factor(group1),font.size=14,title=msig.name,
                                      showCategory=enrich.ncat) + facet_grid(~group2) +
      scale_y_discrete(labels=function(x)
       str_wrap(str_replace_all(str_to_title(tolower(gsub("_"," ",x))),chg.enrich.terms), width=40)) +
      scale_x_discrete(labels=function(x) str_wrap(gsub("_"," ",x),width=10)) +
      scale_color_gradientn(colors=rev(pal(10)),
                            limits=c(0,enrich.pval.co))+
      theme(strip.text=element_text(size=14))
     vals.plot$msig.h
    }) # withProgress
   }) #isolate
  },height=function(){100*react.tab3.enrich.plot.size()})
 
 ## clicking on the export button will generate a pdf file 
 ## containing all plots
 observeEvent( input$export,{
  if(input$export==0) return()
  gen.go <- react.val.gen.go()
  out.loc <- react.tab3.out.loc()
  grp.name <- react.tab3.grp.name()
  withProgress(message="Saving pdf",{
   pad = react.tab3.venn.pad()
   file=file.path(out.loc,"plots.pdf")
   if(file.exists(file)) file.remove(file)
   pdf(file, onefile = TRUE,width=input$pdf.width,height=input$pdf.height)
   if(!is.null(vals.plot$venn.up1))
    gridExtra::grid.arrange(vals.plot$venn.up1,vals.plot$venn.up2,ncol=2,
                           padding=unit(pad,"line"),
                           top=textGrob("Up-regulated Genes",gp=gpar(fontsize=20)),
                           vp=viewport(width=0.8, height=0.8))
   if(!is.null(vals.plot$venn.dwn1))
    gridExtra::grid.arrange(vals.plot$venn.dwn1,vals.plot$venn.dwn2,ncol=2,
                           padding=unit(pad,"line"),
                           top=textGrob("Down-regulated Genes",gp=gpar(fontsize=20)),
                           vp=viewport(width=0.8, height=0.8))
   if(!is.null(vals.plot$hm.up))
    gridExtra::grid.arrange(vals.plot$hm.up,vals.plot$hm.dwn,
                                       ncol=2,padding=unit(pad,"line"),
                           top=textGrob("Significance of overlaps",gp=gpar(fontsize=20)),
                           vp=viewport(width=0.8, height=0.8))
   if(!is.null(vals.plot$upset.up))
    gridExtra::grid.arrange(vals.plot$upset.up,vp=viewport(width=0.8, height=0.8))
   if(!is.null(vals.plot$upset.dwn))
   gridExtra::grid.arrange(vals.plot$upset.dwn,vp=viewport(width=0.8, height=0.8))
   if(!is.null(vals.plot$msig.mf))
    gridExtra::grid.arrange(vals.plot$msig.mf,vp=viewport(width=0.8, height=0.8))
   if(!is.null(vals.plot$msig.bp))
    gridExtra::grid.arrange(vals.plot$msig.bp,vp=viewport(width=0.8, height=0.8))
   if(!is.null(vals.plot$msig.cc))
    gridExtra::grid.arrange(vals.plot$msig.cc,vp=viewport(width=0.8, height=0.8))
   if(!is.null(vals.plot$msig.curate))
    gridExtra::grid.arrange(vals.plot$msig.curate,vp=viewport(width=0.8, height=0.8))
   if(!is.null(vals.plot$msig.h))
    gridExtra::grid.arrange(vals.plot$msig.h,vp=viewport(width=0.8, height=0.8))
   if(!is.null(vals.plot.volcano))
    for(i in 1:length(grp.name))
     gridExtra::grid.arrange(vals.plot.volcano[[paste0("tab3.volcano.grp",i)]],vp=viewport(width=0.8, height=0.8))
    dev.off()
  })#withProgress
 }) # observeEvent
 
 ## clicking on the export button will generate a txt file 
 ## containing filtered gene lists
 observeEvent( input$exportGene,{
  if(input$exportGene==0) return()
  withProgress(message="Exporting gene lists....",{
   grp.name = react.tab3.grp.name()
   grp.up = tab3.grp.up()
   grp.dwn = tab3.grp.dwn()
   fc.co = react.tab3.fc.cutoff()
   fdr.co = react.tab3.fdr.cutoff()
   meanDiff.co = react.tab3.meanDiff.cutoff()
   projdir = react.setup.proj.dir()
   out.loc = react.tab3.out.loc()
   for (i in (1:length(grp.name))){
    # filter the complete table file
    in.file=paste0(gsub("_","",grp.name[i]),".complete.txt")
    df=read.table(paste0(projdir,"outputs/diff_analysis_rslt/tables/",in.file),
                  sep="\t",header=TRUE)
    df.filt=df[grp.up[[i]]$row.pass,]
    df.filt$norm.diff.mean = df[grp.up[[i]]$row.pass,"norm.diff.mean"]
    # write filter up txt
    out.up.file=paste0(gsub("_","",grp.name[i]),
                       paste0("_FC_",fc.co[i],"_FDR_",fdr.co[i],
                              "_meanDiff_",meanDiff.co[i]),".up.txt")
    write.table(df.filt,file=paste0(out.loc,"/",out.up.file),sep="\t",row.names = FALSE,
                quote=FALSE)
    df.filt=df[grp.dwn[[i]]$row.pass,]
    df.filt$norm.diff.mean = df[grp.dwn[[i]]$row.pass,"norm.diff.mean"]
    out.dwn.file=paste0(gsub("_","",grp.name[i]),
                        paste0("_FC_",fc.co[i],"_FDR_",fdr.co[i],
                               "_meanDiff_",meanDiff.co[i]),".down.txt")
    write.table(df.filt,file=paste0(out.loc,"/",out.dwn.file),sep="\t",row.names = FALSE,
                quote=FALSE)
   }
  })#withProgress
 }) # observeEvent
 
 output$tab3.text <- renderPrint({ grp.up = tab3.grp.up(); length(grp.up)})
 
 # Set up the root directory for folder selection
 shinyDirChoose(input, "del_folder", root  = volumes)
 
 # Observe the folder selection and display the selected path
 observe({
  req(input$del_folder)
  
  # Get the selected folder path
  folder_path <- parseDirPath(volumes, input$del_folder)
  
  # Display the folder path
  output$selected_folder <- renderText({
   paste("Selected folder:", folder_path)
  })
  
  # Display merged fastq size
  output$dir_size_merged_fastq <- renderUI({
   HTML(paste("<b>Merged fastq size:", 
         get_dir_size(file.path(folder_path,"outputs","merged_fastq")),"<br>",
         "Trimmed fastq size:",
         get_dir_size(file.path(folder_path,"outputs","trim")),"</b><br><br>"))
  })
  
  # Display bam file size
  output$dir_size_bam <- renderUI({
   HTML(paste("<b>BAM files size:", 
              get_dir_size(file.path(folder_path,"outputs","STAR_2pass")),"</b><br><br>"))
  })
  
 })
 
 # Observe the delete project folder button click
 observeEvent(input$delete, {
  # Show a modal dialog asking for confirmation
  showModal(modalDialog(
   title = "Confirm Deletion",
   "Are you sure you want to delete this folder?",
   easyClose = FALSE,
   footer = tagList(
    modalButton("Cancel"),
    actionButton("confirm_delete", "Yes, Delete")
   )
  ))
 })
 
 # Handle the actual deletion when the user confirms
 observeEvent(input$confirm_delete, {
  # Close the modal dialog
  removeModal()
  req(input$del_folder)

  # Get the selected folder path
  folder_path <- parseDirPath(volumes, input$del_folder)
  
  if (dir.exists(folder_path)) {
   # Attempt to delete the folder
   unlink(folder_path, recursive = TRUE)
   
   if (!dir.exists(folder_path)) {
    output$status <- renderText("Folder deleted successfully.")
   } else {
    output$status <- renderText("Failed to delete the folder.")
   }
  } else {
   output$status <- renderText("Folder does not exist.")
  }
 })
 
 # Observe the delete button click for delete all files except samples.txt
 observeEvent(input$delete_keep_samples, {
  # Show a modal dialog asking for confirmation
  showModal(modalDialog(
   title = "Confirm Deletion",
   "Are you sure you want to delete every files in this folder except samples.txt?",
   easyClose = FALSE,
   footer = tagList(
    modalButton("Cancel"),
    actionButton("confirm_delete_except_samples", "Yes, Delete")
   )
  ))
 })
 
 # Handle the actual deletion when the user confirms
 observeEvent(input$confirm_delete_except_samples, {
  # Close the modal dialog
  removeModal()
  req(input$del_folder)
  
  # Get the selected folder path
  folder_path <- parseDirPath(volumes, input$del_folder)

  if (dir.exists(folder_path)) {
   # Specify the filename you want to keep
   file_to_keep <- "samples.txt"
   
   # Get a list of all files and directories in the specified folder
   all_files <- list.files(folder_path, full.names = TRUE)
   
   # Exclude the file you want to keep
   files_to_delete <- setdiff(all_files, file.path(folder_path, file_to_keep))
   
   print("files_to_delete")
   print(files_to_delete)
   
   # Attempt to delete the folder
   unlink(files_to_delete, recursive = TRUE)
   
   if (length(list.files(folder_path))==1) {
    output$delete_keep_samples_status <- renderText("All files except samples.txt deleted successfully.")
   } else {
    output$delete_keep_samples_status <- renderText("Failed to delete files.")
   }
  } else {
   output$delete_keep_samples_status <- renderText("Folder does not exist.")
  }
 })
 
 # Observe the delete merged fastq files button click
 observeEvent(input$delete_merged_fastq, {
  # Show a modal dialog asking for confirmation
  showModal(modalDialog(
   title = "Confirm Deletion",
   "Are you sure you want to delete merged and/or trimmed fastq files?",
   easyClose = FALSE,
   footer = tagList(
    modalButton("Cancel"),
    actionButton("confirm_delete_merged_fastq", "Yes, Delete")
   )
  ))
 })
 
 # Handle the actual deletion when the user confirms
 observeEvent(input$confirm_delete_merged_fastq, {
  # Close the modal dialog
  removeModal()
  req(input$del_folder)
  
  # Get the selected folder path
  folder_path <- parseDirPath(volumes, input$del_folder)
  merged_fastq_path <- file.path(folder_path,"outputs","merged_fastq")
  trim_fastq_path <- file.path(folder_path,"outputs","trim")
  
  if (dir.exists(merged_fastq_path) || dir.exists(trim_fastq_path)) {
   # Attempt to delete the folder
   unlink(merged_fastq_path, recursive = TRUE)
   # Attempt to delete the folder
   unlink(trim_fastq_path, recursive = TRUE)
   
   if (!dir.exists(merged_fastq_path) || !dir.exists(trim_fastq_path) ) {
    output$delete_merged_fastq_status <- renderText("Merged and/or trimmed fastq files are deleted successfully.")
   } else {
    output$delete_merged_fastq_status <- 
     renderText("Failed to delete merged and/or trimmed fastq files.")
   }
  } else {
   output$delete_merged_fastq_status <- renderText("Merged and/or trimmed fastq files do not exist.")
  }
 })
 
 
 # Observe the delete bam files button click
 observeEvent(input$delete_bam, {
  # Show a modal dialog asking for confirmation
  showModal(modalDialog(
   title = "Confirm Deletion",
   "Are you sure you want to delete BAM files?",
   easyClose = FALSE,
   footer = tagList(
    modalButton("Cancel"),
    actionButton("confirm_delete_bam", "Yes, Delete")
   )
  ))
 })
 
 # Handle the actual deletion when the user confirms
 observeEvent(input$confirm_delete_bam, {
  # Close the modal dialog
  removeModal()
  req(input$del_folder)
  
  # Get the selected folder path
  folder_path <- parseDirPath(volumes, input$del_folder)
  bam_path <- file.path(folder_path,"outputs","STAR_2pass")

  if (dir.exists(bam_path)) {
   # Attempt to delete the folder
   unlink(bam_path, recursive = TRUE)

   if (!dir.exists(bam_path)) {
    output$delete_bam_status <- renderText("BAM files are deleted successfully.")
   } else {
    output$delete_bam_status <- 
     renderText("Failed to BAM files.")
   }
  } else {
   output$delete_bam_status <- renderText("BAM files do not exist.")
  }
 })
}

# app
shinyApp(ui = ui, server = server)
