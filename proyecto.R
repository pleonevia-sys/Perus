# Cargamos las librerias
library(MASS)
library(ggplot2)
library(multimode)
library(Rcpp)
library(pkgbuild)



setwd("~/unam/Simulacion_Estocastica")
sourceCpp("Muestreo_Gibbs.cpp")

# Cargamos el dataset de galaxias
data(galaxies)
data(acidity)
data(diamonds)

# Asignamos a un vector
datos1 <- as.vector(galaxies)
datos2 <- as.vector(acidity)
datos3 <- diamonds$carat

# Ver un resumen
print(summary(datos1))
print(summary(datos2))
print(summary(datos_diamantes))

# Graficar el histograma inicial para visualizar las múltiples modas (clústeres)
ggplot(data.frame(x = datos1), aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 60, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkred", size = 1) +
  labs(title = "Distribución de Velocidades de Galaxias",
       x = "Velocidad (10^3 km/s)",
       y = "Densidad") +
  theme_minimal()




ggplot(data.frame(x = datos2), aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 60, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkred", size = 1) +
  labs(title = "Distribución de Velocidades de Galaxias",
       x = "Acidez",
       y = "Densidad") +
  theme_minimal()



ggplot(data.frame(x = datos3), aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 60, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkred", size = 1) +
  labs(title = "Distribución de Velocidades de Galaxias",
       x = "Acidez",
       y = "Densidad") +
  theme_minimal()



# Función genérica para inicializar el modelo DPGMM
inicializar_modelo_dpgmm <- function(datos, K_max = 30) {
  n <- length(datos)
  
  # 1. Asignaciones iniciales aleatorias z_i para cada dato (de 1 a K_max)
  z <- sample(1:K_max, n, replace = TRUE)
  
  # 2. Inicializar medias de los clústeres basadas en la distribución de los datos
  mu <- rnorm(K_max, mean = datos, sd = sd(datos) / 2)
  
  # 3. Inicializar varianzas para cada clúster (ej. Inverse-Gamma o varianza común inicial)
  sigma2 <- rep(var(datos) / 2, K_max)
  
  # 4. Inicializar variables Beta auxiliares para el Stick-Breaking
  # Usamos un valor inicial para alpha
  alpha_val <- 1.0
  V <- rbeta(K_max - 1, 1, alpha_val)
  V[K_max] <- 1 # Asegurar que el último sume 1
  
  # Calcular los pesos pi a partir de V
  pi_pesos <- numeric(K_max)
  pi_pesos[1] <- V[1]
  for (k in 2:K_max) {
    pi_pesos[k] <- V[k] * prod(1 - V[1:(k-1)])
  }
  pi_pesos <- pi_pesos / sum(pi_pesos) # Normalizar por seguridad numérica
  
  # Retornar la lista con el estado inicial del modelo
  list(
    z = z,
    mu = mu,
    sigma2 = sigma2,
    pi = pi_pesos,
    V = V,
    alpha = alpha_val,
    K_max = K_max,
    n = n
  )
}





# Función para ejecutar el Gibbs Sampler completo de manera genérica
ejecutar_gibbs_dpgmm <- function(datos, n_iter = 1000, K_max = 30) {
  
  # 1. Inicializar el modelo
  modelo <- inicializar_modelo_dpgmm(datos, K_max)
  n <- modelo$n
  
  # Hiperparámetros previos (priors) por defecto ajustados a los datos
  m0 <- mean(datos)
  tau0_sq <- var(datos)
  a_sigma <- 2.0
  b_sigma <- var(datos)
  
  # Hiperparámetros para alpha ~ Gamma(a_alpha, b_alpha)
  a_alpha <- 2.0
  b_alpha <- 2.0
  
  # Matrices para almacenar los resultados (muestras de la cadena MCMC)
  historial_mu <- matrix(NA, nrow = n_iter, ncol = K_max)
  historial_alpha <- numeric(n_iter)
  historial_k_activos <- numeric(n_iter)
  
  # 2. Bucle principal de MCMC
  for (iter in 1:n_iter) {
    
    # --- PASO 1 (OPTIMIZADO EN C++): Actualizar asignaciones z_i ---
    # Aquí llamamos directamente a la función ultrarrápida compilada desde C++
    modelo$z <- actualizar_z_cpp(datos, modelo$mu, modelo$sigma2, modelo$pi, K_max)
    
    # --- PASO 2 y 3: Actualizar Medias (mu) y Varianzas (sigma2) por clúster ---
    for (k in 1:K_max) {
      puntos_k <- datos[modelo$z == k]
      n_k <- length(puntos_k)
      
      if (n_k > 0) {
        media_datos <- mean(puntos_k)
        
        # Posteriores para la media mu_k
        v_post <- 1 / ( (1 / tau0_sq) + (n_k / modelo$sigma2[k]) )
        m_post <- v_post * ( (m0 / tau0_sq) + (n_k * media_datos / modelo$sigma2[k]) )
        modelo$mu[k] <- rnorm(1, mean = m_post, sd = sqrt(v_post))
        
        # Posteriores para la varianza sigma2_k (Inverse-Gamma)
        a_post <- a_sigma + n_k / 2
        b_post <- b_sigma + sum((puntos_k - modelo$mu[k])^2) / 2
        modelo$sigma2[k] <- 1 / rgamma(1, shape = a_post, rate = b_post)
        
      } else {
        modelo$mu[k] <- rnorm(1, mean = m0, sd = sqrt(tau0_sq))
        modelo$sigma2[k] <- 1 / rgamma(1, shape = a_sigma, rate = b_sigma)
      }
    }
    
    # --- PASO 4: Actualizar Pesos (Stick-Breaking) ---
    # Contar cuántos datos cayeron en cada clúster
    n_k_counts <- tabulate(modelo$z, nbins = K_max)
    
    # Muestrear las variables Beta auxiliares V
    for (k in 1:(K_max - 1)) {
      suma_posteriores <- sum(n_k_counts[(k + 1):K_max])
      modelo$V[k] <- rbeta(1, shape1 = 1 + n_k_counts[k], shape2 = modelo$alpha + suma_posteriores)
    }
    modelo$V[K_max] <- 1
    
    # Reconstruir los pesos pi a partir de V
    modelo$pi[1] <- modelo$V[1]
    for (k in 2:K_max) {
      modelo$pi[k] <- modelo$V[k] * prod(1 - modelo$V[1:(k-1)])
    }
    modelo$pi <- modelo$pi / sum(modelo$pi)
    
    # --- PASO 5: Actualizar el Parámetro de Concentración desconocido (alpha) ---
    # Usando el método auxiliar de Escobar & West (1995)
    k_activos <- length(unique(modelo$z))
    eta <- rbeta(1, shape1 = modelo$alpha + 1, shape2 = n)
    pi_eta <- (a_alpha + k_activos - 1) / ( (a_alpha + k_activos - 1) + n * (b_alpha - log(eta)) )
    
    if (runif(1) < pi_eta) {
      modelo$alpha <- rgamma(1, shape = a_alpha + k_activos, rate = b_alpha - log(eta))
    } else {
      modelo$alpha <- rgamma(1, shape = a_alpha + k_activos - 1, rate = b_alpha - log(eta))
    }
    
    # Guardar trazas para diagnósticos
    historial_mu[iter, ] <- modelo$mu
    historial_alpha[iter] <- modelo$alpha
    historial_k_activos[iter] <- k_activos
    
    if (iter %% 500 == 0) {
      cat(sprintf("Iteración %d/%d completada. Clústeres activos: %d | Alpha: %.2f\n", 
                  iter, n_iter, k_activos, modelo$alpha))
    }
  }
  
  # Retornar historial y estado final
  list(
    estado_final = modelo,
    historial_mu = historial_mu,
    historial_alpha = historial_alpha,
    historial_k_activos = historial_k_activos
  )
}




# Prueba con el dataset de galaxias
resultado_galaxias <- ejecutar_gibbs_dpgmm(datos = datos1)
resultado_acidez <- ejecutar_gibbs_dpgmm(datos = datos2)
resultado_diamante <- ejecutar_gibbs_dpgmm(datos = datos3)


# 1. Diagnósticos de convergencia (¿Se estabilizó la cadena?)
par(mfrow = c(1, 2))

# Gráfica de la evolución del número de clústeres activos
plot(resultado_galaxias$historial_k_activos, type = "l", col = "darkblue",
     main = "Evolución de Clústeres Activos", xlab = "Iteración", ylab = "Número de Clústeres (k*)")

# Gráfica de la evolución del parámetro de concentración (alpha)
plot(resultado_galaxias$historial_alpha, type = "l", col = "darkred",
     main = "Evolución de Alpha", xlab = "Iteración", ylab = "Valor de Alpha")

par(mfrow = c(1, 1)) # Regresar a una sola ventana gráfica


# 2. Estimación de la Densidad Posterior
# Descartamos el 20% inicial de la cadena como "burn-in" (calentamiento)
n_iter_total <- length(resultado_galaxias$historial_alpha)
burn_in <- floor(n_iter_total * 0.2)
indices_post <- (burn_in + 1):n_iter_total

# Generar una malla de puntos para evaluar la curva de densidad estimada
x_grid <- seq(min(datos1) - 5, max(datos1) + 5, length.out = 200)
densidad_estimada <- numeric(length(x_grid))

# Promediar las densidades gaussianas de las muestras posteriores seleccionadas
for (t in indices_post) {
  # Extraer medias, varianzas y pesos de la iteración t del historial
  mu_t <- resultado_galaxias$historial_mu[t, ]
  # Nota: si guardaste sigma2 en el historial puedes usarlo; aquí usamos una aproximación basada en el estado final
  sigma_t <- sqrt(resultado_galaxias$estado_final$sigma2) 
  pi_t <- resultado_galaxias$estado_final$pi
  
  # Calcular la mezcla de gaussianas para esta iteración
  dens_t <- rep(0, length(x_grid))
  for (k in 1:length(mu_t)) {
    if (!is.na(mu_t[k]) && pi_t[k] > 1e-4) {
      dens_t <- dens_t + pi_t[k] * dnorm(x_grid, mean = mu_t[k], sd = sigma_t[k])
    }
  }
  densidad_estimada <- densidad_estimada + dens_t
}

# Normalizar dividiendo entre el número de iteraciones usadas
densidad_estimada <- densidad_estimada / length(indices_post)

# 3. Graficar el resultado final: Datos reales vs. Curva de Densidad Estimada por el DPGMM
hist(datos1, probability = TRUE, breaks = 15, col = "gray90", border = "white",
     main = "Estimación de Densidad con DPMM (Gibbs Sampler)",
     xlab = "Datos", ylab = "Densidad")
lines(x_grid, densidad_estimada, col = "darkred", lwd = 3, lty = 2)









# 2. Estimación de la Densidad Posterior
# Descartamos el 20% inicial de la cadena como "burn-in" (calentamiento)
n_iter_total <- length(resultado_acidez$historial_alpha)
burn_in <- floor(n_iter_total * 0.2)
indices_post <- (burn_in + 1):n_iter_total

# Generar una malla de puntos para evaluar la curva de densidad estimada
x_grid <- seq(min(datos2) - 5, max(datos2) + 5, length.out = 200)
densidad_estimada <- numeric(length(x_grid))

# Promediar las densidades gaussianas de las muestras posteriores seleccionadas
for (t in indices_post) {
  # Extraer medias, varianzas y pesos de la iteración t del historial
  mu_t <- resultado_acidez$historial_mu[t, ]
  # Nota: si guardaste sigma2 en el historial puedes usarlo; aquí usamos una aproximación basada en el estado final
  sigma_t <- sqrt(resultado_acidez$estado_final$sigma2) 
  pi_t <- resultado_acidez$estado_final$pi
  
  # Calcular la mezcla de gaussianas para esta iteración
  dens_t <- rep(0, length(x_grid))
  for (k in 1:length(mu_t)) {
    if (!is.na(mu_t[k]) && pi_t[k] > 1e-4) {
      dens_t <- dens_t + pi_t[k] * dnorm(x_grid, mean = mu_t[k], sd = sigma_t[k])
    }
  }
  densidad_estimada <- densidad_estimada + dens_t
}

# Normalizar dividiendo entre el número de iteraciones usadas
densidad_estimada <- densidad_estimada / length(indices_post)

# 3. Graficar el resultado final: Datos reales vs. Curva de Densidad Estimada por el DPGMM
hist(datos2, probability = TRUE, breaks = 15, col = "gray90", border = "white",
     main = "Estimación de Densidad con DPMM (Gibbs Sampler)",
     xlab = "Datos", ylab = "Densidad")
lines(x_grid, densidad_estimada, col = "darkred", lwd = 3, lty = 2)






# 2. Estimación de la Densidad Posterior
# Descartamos el 20% inicial de la cadena como "burn-in" (calentamiento)
n_iter_total <- length(resultado_diamante$historial_alpha)
burn_in <- floor(n_iter_total * 0.2)
indices_post <- (burn_in + 1):n_iter_total

# Generar una malla de puntos para evaluar la curva de densidad estimada
x_grid <- seq(min(datos3) - 5, max(datos3) + 5, length.out = 200)
densidad_estimada <- numeric(length(x_grid))

# Promediar las densidades gaussianas de las muestras posteriores seleccionadas
for (t in indices_post) {
  # Extraer medias, varianzas y pesos de la iteración t del historial
  mu_t <- resultado_diamante$historial_mu[t, ]
  # Nota: si guardaste sigma2 en el historial puedes usarlo; aquí usamos una aproximación basada en el estado final
  sigma_t <- sqrt(resultado_diamante$estado_final$sigma2) 
  pi_t <- resultado_diamante$estado_final$pi
  
  # Calcular la mezcla de gaussianas para esta iteración
  dens_t <- rep(0, length(x_grid))
  for (k in 1:length(mu_t)) {
    if (!is.na(mu_t[k]) && pi_t[k] > 1e-4) {
      dens_t <- dens_t + pi_t[k] * dnorm(x_grid, mean = mu_t[k], sd = sigma_t[k])
    }
  }
  densidad_estimada <- densidad_estimada + dens_t
}

# Normalizar dividiendo entre el número de iteraciones usadas
densidad_estimada <- densidad_estimada / length(indices_post)

# 3. Graficar el resultado final: Datos reales vs. Curva de Densidad Estimada por el DPGMM
hist(datos3, probability = TRUE, breaks = 60,
     main = "Estimación de Densidad con DPMM (Gibbs Sampler)",
     xlab = "Datos", ylab = "Densidad")
lines(x_grid, densidad_estimada, col = "darkred", lwd = 3, lty = 2)

