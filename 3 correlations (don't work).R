
test <- df # type
test$EntryNo <- NULL
test$TaxonBioStatus <- NULL
test$TaxonGrowthForm <- NULL
test$`Overhanging Ring` <- NULL
test$Perc <- NULL
test$Weight <- NULL
test$GroundCover <- NULL



test.wide <- test %>%
  pivot_wider(
    id_cols = c(Year, Plot),
    names_from = `Verbatim Species`,
    values_from = Proportion,
    values_fn = mean, # if there are multiple values average them
    values_fill = 0
  )

test.wide <- as.data.frame(test.wide)
test.wide[is.na(test.wide)] <- 0 

species.cor <- test.wide %>%
  group_by(Year, Plot) 

species.cor$Year <- NULL
species.cor$Plot <- NULL

species.cor <- species.cor %>%
  select(where(~ sd(.x, na.rm = TRUE) > 0)) %>%
  cor(use = "pairwise.complete.obs") %>% 
  as.data.frame()


species.cor <- species.cor %>% rownames_to_column()
species.cor <- species.cor %>% rename(Species = rowname)


brapin <- species.cor[, c("Species", "BRAPIN")]

brapin <- brapin %>% arrange(BRAPIN)
brapin <- brapin[1:(nrow(brapin)-1), ]
brapin$direction <- ifelse(brapin$BRAPIN <0, "Neg", "Pos")

ggplot()+
  theme_bw()+
  geom_col(data = brapin, aes(x= reorder(Species, BRAPIN) , 
                              y = BRAPIN, fill = direction))+
  scale_fill_manual(values = c("purple", "red"))+
  theme(axis.text = element_text(angle = 90, hjust =1, size =6))


# all scientific names
sort(unique(test$NVSSpeciesName))

# Individual correlations
names(df)

my.selection <- df %>%  filter(NVSSpeciesName %in% 
                                 c("Viola filicaulis", 
                                   x = "Brachyscome pinnata"))
my.selection.wide <- my.selection %>%
  select(Year, Plot, NVSSpeciesName, Proportion) %>%
  group_by(Year, Plot)%>%
  pivot_wider(names_from =  NVSSpeciesName,
              values_from = Proportion,
              values_fill = 0,
              values_fn = mean)  # duplicate settings

head(my.selection.wide)

ggplot()+
  geom_point(data = my.selection.wide, aes(y = `Viola filicaulis`, 
                                           x = `Brachyscome pinnata`))

