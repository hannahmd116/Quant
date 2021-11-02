######## Get necessary packages ##########
{
if (!require(behavr)) {install.packages('behavr');library(behavr)}
if (!require(damr, quietly=T)) {install.packages('damr');library(damr)}
if (!require(ggetho, quietly=T)) {install.packages('ggetho');library(ggetho)}
if (!require(sleepr, quietly=T)) {install.packages('sleepr');library(sleepr)}
if (!require(zeitgebr, quietly=T)) {install.packages('zeitgebr');library(zeitgebr)}
if (!require(svDialogs, quietly=T)) {install.packages('svDialogs');require(svDialogs)}
if (!require(dplyr, quietly=T)) {install.packages('dplyr');require(dplyr)}
}

######### Choose experimental folder ########
{
wd <- dlg_dir(title="Choose the experiment folder")$res
setwd(wd)
}

######## Read in data and set up aesthetics #########
### read in metadata and link
{
metafile <- dlg_open(default = wd, title="Choose the metadata file")$res
metadata <- fread(metafile)
metadata <- link_dam_metadata(metadata, result_dir = wd)
### set up aesthetics
by_colour <- dlg_input(message = "Which variable do you want to colour by? CASE SENSITIVE\n(For e.g. genotype, driver, age, treatment, etc.)", 
                        default="genotype")$res
by_facet <- dlg_input(message = "Which variable do you want to create multiple panels for? CASE SENSITIVE\n(For e.g. genotype, driver, age, treatment, etc.)", 
                      default="")$res
}

######## Load DAM data and choose what to include #########
{
### filter what to be included
by_colour_list <- metadata %>% select(contains(by_colour)) %>% unique() %>% pull()
if (length(by_facet)==0){
  to_plot <- dlg_list(choices = by_colour_list, multiple=T, title="Which one do you want to include to analyse?")$res
  metadata2 <- metadata[(get(by_colour) %in% to_plot)]
} else {
  by_facet_list <- metadata %>% select(contains(by_facet)) %>% unique() %>% pull()
  list <- c(by_colour_list, by_facet_list)
  to_plot <- dlg_list(choices = list, multiple=T, title="Which one do you want to include to analyse?")$res
  metadata2 <- metadata[(get(by_colour) %in% to_plot)&(get(by_facet) %in% to_plot)]
}
dt <- load_dam(metadata2, FUN = sleepr::sleep_dam_annotation)
summary(dt)
ggetho(dt, aes(z=asleep)) +
  stat_ld_annotations(height = 1) +
  stat_tile_etho() +
  labs(fill="sleep value")
}

####### (Optional) Automatic curation of dead animals #######
{
dt_original <- dt
dt <- curate_dead_animals(dt)
dead_flies <- setdiff(dt_original[, meta=T],dt[, meta=T]) %>% select(-file_info)
fwrite(dead_flies, "dead_flies.csv")
summary(dt)
ggetho(dt, aes(z=asleep)) +
  stat_ld_annotations(ypos = "top")+
  stat_tile_etho() +
  labs(fill="sleep value")
}

##### (Optional) Filter by how many days they survived #####
{
day_cutoff <- as.numeric(dlg_input(message="At least how many days do the animals need to be alive to be included?", default="1")$res)
lifespan_dt <- dt[, .(lifespan = max(t)), by=id]        # calculate lifespan for each animals
valid_ids <- lifespan_dt[lifespan >= days(day_cutoff), id]                # filter for desirable lifespan
decision <- dlg_input(message=paste0(length(valid_ids)," animals meet the criteria, do you want to continue?\nY or N?"), 
                        default="Y")$res
  if (decision == "Y") {
    dt <- dt[id %in% valid_ids]    # apply this filter
  } else { dt <- dt}
  
  summary(dt)
  ggetho(dt, aes(z=asleep)) +
    stat_ld_annotations(ypos = "top")+
    stat_tile_etho() +
    labs(fill="sleep value")
}
 
###### if we want to look at only certain time period
dt <- dt[t %between% c(days(0), days(2))]
ggetho(dt, aes(z=asleep)) +
  stat_ld_annotations(ypos = "top")+
  stat_tile_etho() +
  labs(fill="sleep value")

 
####### Sleep traces ##### 
{
  no_of_days <- max(dt$t)/hours(24)
  width <- no_of_days*6
  sleep_line_days <- ggetho(dt, aes_string(y="asleep", colour=by_colour)) +
     stat_pop_etho() +
     stat_ld_annotations() +
     facet_grid(c(by_facet)) +
     scale_y_continuous(name= "Fraction of time sleeping",labels = scales::percent) +
     theme_bw() +
     theme(legend.position="bottom")
   ggsave("Sleep_line_days.pdf", sleep_line_days, width=8.1, height=6.11)
   
   sleep_line <- ggetho(dt, aes_string(y="asleep", colour=by_colour), time_wrap = hours(24)) +
     stat_pop_etho() +
     stat_ld_annotations() +
     facet_grid(c(by_facet)) +
     scale_y_continuous(name= "Fraction of time sleeping",labels = scales::percent) +
     theme_bw() +
     theme(legend.position="bottom")
   ggsave("Sleep_line_avergae.pdf", sleep_line, width=8, height=6.11)
 }
 
####### Sleep analysis with phase #####
{
dt_phase <- dt[, phase := ifelse(t %% hours(24) < hours(12), "L", "D")]
summary_dt_phase <- rejoin(dt_phase[,.(fraction_all = mean(asleep),
                                       fraction_light = mean(asleep[phase == "L"]),
                                       fraction_dark = mean(asleep[phase == "D"]),
                                       amount_all = sum(asleep)/no_of_days,
                                       amount_light = sum(asleep[phase == "L"])/no_of_days,
                                       amount_dark = sum(asleep[phase == "D"])/no_of_days
                                      ),by=id])
no_col <- ncol(summary_dt_phase)
summary_dt_melted <- reshape(summary_dt_phase, varying=(no_col-5):no_col, sep="_", direction="long", timevar="phase")

fraction_phase <- ggplot(summary_dt_melted, aes_string(x="phase", y="fraction", fill=by_colour)) +
  geom_boxplot(outlier.colour = NA, width=0.5, size=0.1, position=position_dodge(w=0.5)) +
  geom_jitter(alpha=.6, size=0.6, position=position_jitterdodge(dodge.width=0.5, jitter.width=0.05)) +
  facet_grid(c(by_facet)) +
  scale_y_continuous(name= "Fraction of time sleeping",labels = scales::percent) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("Sleep_fraction_phase.pdf",fraction_phase)

amount_phase <- ggplot(summary_dt_melted, aes_string(x="phase", y="amount", fill=by_colour)) +
  geom_boxplot(outlier.colour = NA, width=0.5, size=0.1, position=position_dodge(w=0.5)) +
  geom_jitter(alpha=.6, size=0.6, position=position_jitterdodge(dodge.width=0.5, jitter.width=0.05)) +
  facet_grid(c(by_facet)) +
  scale_y_continuous(name= "Amount of sleeping (minutes)") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank())
ggsave("Sleep_amount_phase.pdf",amount_phase)

### Bout analysis
bout_dt <- bout_analysis(asleep, dt)
bout_dt <- bout_dt[asleep == TRUE, -"asleep"]
bout_dt_phase <- bout_dt[, phase := ifelse(t %% hours(24) < hours(12), "L", "D")]
bout_summary_phase <- rejoin(bout_dt_phase[,.(bouts_all = length(duration)/no_of_days,length_all = mean(duration)/(60*no_of_days),
                                              bouts_light = length(duration[phase == "L"])/no_of_days,
                                              length_light = mean(duration[phase == "L"])/(60*no_of_days),
                                              bouts_dark = length(duration[phase == "D"])/no_of_days,
                                              length_dark = mean(duration[phase == "D"])/(60*no_of_days)
                                              ),by=id])
no_col_bout <- ncol(bout_summary_phase)
summary_bout_melted <- reshape(bout_summary_phase, varying=(no_col_bout-5):no_col_bout, sep="_", direction="long", timevar="phase")

bout_number <- ggplot(summary_bout_melted, aes_string(x="phase", y="bouts", fill=by_colour)) + 
  geom_boxplot(outlier.colour = NA, width=0.5, size=0.1, position=position_dodge(w=0.5)) +
  geom_jitter(alpha=.6, size=0.6, position=position_jitterdodge(dodge.width=0.5, jitter.width=0.05)) +
  facet_grid(c(by_facet)) +
  scale_y_continuous(name= "Number of sleep bouts") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("Bouts_number_phase.pdf", bout_number)

bout_length <- ggplot(summary_bout_melted, aes_string(x="phase", y="length", fill=by_colour)) + 
  geom_boxplot(outlier.colour = NA, width=0.5, size=0.1, position=position_dodge(w=0.5)) +
  geom_jitter(alpha=.6, size=0.6, position=position_jitterdodge(dodge.width=0.5, jitter.width=0.05)) +
  facet_grid(c(by_facet)) +
  scale_y_continuous(name= "Average length of sleep bouts") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("Bouts_length_phase.pdf", bout_length)
}

####### Output data as a csv file #####
{
overall_summary <- merge.data.frame(summary_dt_melted, summary_bout_melted, all=T) %>%
  select(-file_info,-region_id,-experiment_id,-start_datetime,-stop_datetime)
fwrite(overall_summary, "Summary_sleep.csv")
}
