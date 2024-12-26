# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tmap)
library(sf)
library(lmtest)
library(spdep)
library(spatialreg)
library(geosphere)
library(car)

# Membaca data
df <- read_excel("/Users/Golda/College/sem 5/spasial/Tugas 4/data regresi spasial.xlsx")
head(df)

# Membaca shapefile
shapefile <- st_read("/Users/Golda/College/sem 5/spasial/Tugas 4/shapefile jabar/Jawa_Barat_ADMIN_BPS.shp")
shapefile <- shapefile[shapefile$Kabupaten != "Waduk Cirata", ]
colnames(df)[colnames(df) == "Kabupaten/Kota"] <- "Kabupaten"
shapefile <- shapefile %>%left_join(df, by = "Kabupaten")

# Statistika deskriptif
summary(df)

# Normalisasi data untuk menghindari perbedaan skala
df <- df %>%
  mutate(across(c(PDRB, APS, TPAK, KP, ASL, TK, JAK), ~ scale(.)))

# Visualisasi data
variables <- c("TPT", "PDRB", "APS", "TPAK", "KP", "ASL", "TK", "JAK")

# Iterasi untuk setiap variabel
for (var in variables) {
  # Peta
  p <- ggplot(data = shapefile) +
    geom_sf(aes_string(fill = var), color = "black", size = 0.2) +
    scale_fill_gradient(low = "yellow", high = "dark green", name = var) +
    labs(title = paste("Peta Nilai", var, "Per Wilayah")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", vjust = 2, size = 20), # Posisi dan gaya judul
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10) # Margin judul
    )
  print(p) # Tampilkan plot
  
  # Boxplot
  boxplot(df[[var]],
          main = paste("Boxplot", var),
          ylab = var,
          col = "lightblue")
}


# Model OLS untuk regresi linear
ols_model <- lm(TPT ~ PDRB + APS + TPAK + KP + ASL + TK + JAK, data = df)
summary(ols_model)

# Uji asumsi model OLS
  # 1. Normalitas residual
  qqnorm(residuals(ols_model))
  qqline(residuals(ols_model), col = "red")
  shapiro_test <- shapiro.test(residuals(ols_model))
  cat("Uji Shapiro-Wilk:", shapiro_test$p.value, "\n")
  
  # 2. Heteroskedastisitas
  bptest_ols <- bptest(ols_model)
  cat("Uji Breusch-Pagan:", bptest_ols$p.value, "\n")
  
  # 3. Autokorelasi residual
  dw_test <- dwtest(ols_model)
  cat("Uji Durbin-Watson:", dw_test$p.value, "\n")
  
  # 4. Multikolinearitas
  vif_values <- vif(ols_model)
  print(vif_values)

  # 5. Plot residual spasial model terbaik 
  residuals_best <- residuals(ols_model) 
  shapefile$Residuals <- residuals_best
  tm_shape(shapefile) +
    tm_polygons("Residuals", palette = "RdYlBu", midpoint = 0, title = paste("Residuals:", ols_model)) +
    tm_layout(
      outer.margins = c(0.1, 0, 0, 0), # Menambahkan margin atas
      main.title = paste("Residuals Map for OLS"), # Judul di luar kotak
      main.title.size = 1.5, # Ukuran huruf judul
      main.title.fontface = "bold", # Huruf tebal
      main.title.position = "center" # Posisi di tengah
    )  

# Membuat matriks koordinat
coords <- cbind(df$Longitude, df$Latitude)

# Menghitung jarak antar koordinat menggunakan geosphere
distances <- distm(coords)

# Menghitung threshold jarak untuk hubungan spasial
threshold <- quantile(distances[distances > 0], 0.2)  # Naikkan jika terlalu jarang

# Membuat matriks bobot spasial berdasarkan threshold
weights_matrix <- distances < threshold
diag(weights_matrix) <- 0  # Tidak ada hubungan dengan dirinya sendiri
listw <- mat2listw(weights_matrix, style = "W")

# Uji Moran's I 
  # 1. Untuk variabel dependen
  moran_test_dep <- moran.test(df$TPT, listw)
  print(moran_test_dep)
  
  # 2. Untuk variabel independen
  indep_vars <- c("PDRB", "APS", "TPAK", "KP", "ASL", "TK", "JAK")
  
  # Iterasi untuk setiap variabel independen
  for (var in indep_vars) {
    cat("\nUji Moran's I untuk", var, ":\n")
    print(moran.test(df[[var]], listw))
  }
  
  # 3. Untuk residual OLS
  moran_test_resid <- moran.test(residuals(ols_model), listw)
  print(moran_test_resid)

# Model Regresi Spasial
models <- list(
  SLX = lmSLX(TPT ~ PDRB + APS + TPAK + KP + ASL + TK + JAK, data = df, listw = listw),
  SLM = lagsarlm(TPT ~ PDRB + APS + TPAK + KP + ASL + TK + JAK, data = df, listw = listw, zero.policy = TRUE),
  SEM = errorsarlm(TPT ~ PDRB + APS + TPAK + KP + ASL + TK + JAK, data = df, listw = listw, zero.policy = TRUE),
  GSM = sacsarlm(TPT ~ PDRB + APS + TPAK + KP + ASL + TK + JAK, data = df, listw = listw, zero.policy = TRUE),
  SDM = lagsarlm(TPT ~ PDRB + APS + TPAK + KP + ASL + TK + JAK, data = df, listw = listw, type = "mixed", zero.policy = TRUE),
  SDEM = errorsarlm(TPT ~ PDRB + APS + TPAK + KP + ASL + TK + JAK, data = df, listw = listw, etype = "emixed", zero.policy = TRUE),
  GNSM = sacsarlm(TPT ~ PDRB + APS + TPAK + KP + ASL + TK + JAK, data = df, listw = listw, type = "general", zero.policy = TRUE)
)

# Iterasi untuk setiap model
for (model_name in names(models)) {
  model <- models[[model_name]]
  
  # Tampilkan ringkasan model
  cat("\n=== Model:", model_name, "===\n")
  print(summary(model))
  
  # Residual model
  residuals_model <- residuals(model)
  
  # KS-test
  ks_test <- ks.test(residuals_model, "pnorm", mean(residuals_model), sd(residuals_model))
  cat("\nHasil KS-test untuk", model_name, ":\n")
  print(ks_test)
  
  # BP-test
  bp_test <- bptest(residuals_model ~ PDRB + APS + TPAK + KP + ASL + TK + JAK, data = df)
  cat("\nHasil BP-test untuk", model_name, ":\n")
  print(bp_test)
  
  # Moran's I
  moran_test <- moran.test(residuals_model, listw)
  cat("\nHasil Moran's I untuk", model_name, ":\n")
  print(moran_test)
  
  # Menambahkan residuals ke shapefile
  shapefile$Residuals <- residuals_model
  
  # Plot residual spasial pada peta
  cat("\nVisualisasi residual spasial untuk", model_name, ":\n")
  print(
    tm_shape(shapefile) +
      tm_polygons("Residuals", palette = "RdYlBu", midpoint = 0, title = paste("Residuals:", model_name)) +
      tm_layout(
        outer.margins = c(0.1, 0, 0, 0), # Menambahkan margin atas
        main.title = paste("Residuals Map for", model_name), # Judul di luar kotak
        main.title.size = 1.5, # Ukuran huruf judul
        main.title.fontface = "bold", # Huruf tebal
        main.title.position = "center" # Posisi di tengah
      )
  )
}


# Filter hanya model SLM dan SDM dari list models
filtered_models <- models[c("SLM", "SDM")]

# (berdasarkan uji Moran's I: terdapat autokorelasi pada variabel dependen dan independen, sehingga model yang digunakan adalah SLM dan SDM)
# Membandingkan nilai AIC untuk SLM dan SDM
aic_values <- sapply(filtered_models, AIC)
cat("Nilai AIC untuk model SLM dan SDM:\n")
print(aic_values)

# Menentukan model terbaik berdasarkan AIC
best_model <- names(aic_values)[which.min(aic_values)]
cat("\nModel terbaik berdasarkan AIC adalah:", best_model, "\n")

# Pada model SDM terdpat variabel yang tidak signifikan, yaitu APS dan JAK. Variabel APS akan diabaikan  
sdm2_model <- lagsarlm(TPT ~ PDRB + TPAK + KP + ASL + TK + JAK, data = df, listw = listw, type = "mixed", zero.policy = TRUE) 
summary(sdm2_model)

  # KS-test
  ks_test <- ks.test(residuals(sdm2_model), "pnorm", mean(residuals(sdm2_model)), sd(residuals(sdm2_model)))
  cat("Hasil KS-test:\n")
  print(ks_test)
  
  # BP-test 
  bp_test <- bptest(residuals(sdm2_model) ~ PDRB + APS + TPAK + KP + ASL + TK + JAK, data = df)
  cat("\nHasil BP-test:\n")
  print(bp_test)
  
  # Moran's I
  moran_test <- moran.test(residuals(sdm2_model), listw)
  cat("\nHasil Moran's I:\n")
  print(moran_test)

# Plot residual spasial model terbaik 
residuals_best <- residuals(sdm2_model) 

# Menambahkan residuals ke shapefile 
shapefile$Residuals <- residuals_best
tm_shape(shapefile) +
  tm_polygons("Residuals", palette = "RdYlBu", midpoint = 0, title = paste("Residuals:", sdm2_model)) +
  tm_layout(
    outer.margins = c(0.1, 0, 0, 0), # Menambahkan margin atas
    main.title = paste("Residuals Map for SDM"), # Judul di luar kotak
    main.title.size = 1.5, # Ukuran huruf judul
    main.title.fontface = "bold", # Huruf tebal
    main.title.position = "center" # Posisi di tengah
  )

# Menghitung efek langsung, tidak langsung, dan total untuk SDM
impacts_sdm <- impacts(sdm_model, listw = listw, R = 1000)  # Bootstrap untuk interval kepercayaan
summary(impacts_sdm)

# Ekstrak data dari impacts_sdm$res
effects_direct <- impacts_sdm$res$direct
effects_indirect <- impacts_sdm$res$indirect
effects_total <- impacts_sdm$res$total

# Gabungkan efek menjadi matriks
effects_matrix <- cbind(Direct = effects_direct, Indirect = effects_indirect, Total = effects_total)

# Pastikan matriks memiliki nama baris (variabel independen)
rownames(effects_matrix) <- attr(impacts_sdm, "bnames")  # Nama variabel independen

# Visualisasi efek langsung, tidak langsung, dan total
barplot(
  t(effects_matrix),  # Transpose matriks agar variabel menjadi sumbu x
  beside = TRUE,
  legend = rownames(t(effects_matrix)),  # Tambahkan legenda
  col = c("blue", "green", "red"),  # Warna batang
  names.arg = rownames(effects_matrix),  # Nama variabel independen
  main = "Efek Langsung, Tidak Langsung, dan Total",
  xlab = "Variabel Independen",
  ylab = "Magnitudo Efek"
)