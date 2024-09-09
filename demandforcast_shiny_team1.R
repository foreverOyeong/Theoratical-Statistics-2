

library(fresh)
library(DT)

library(dashboardthemes)
library(shiny)
library(shinydashboard)
library(ggplot2)
library(tidyverse)
library(lubridate)

### Themes ##############
mytheme <- create_theme(
  adminlte_color(
    light_blue = '#4CAF50'
  ),
  adminlte_sidebar(
    dark_bg = '#1B5E20',
    dark_hover_bg = '#66BB6A',
    dark_color = '#FFFFFF'
  ),
  adminlte_global(
    content_bg = '#FFFFFF',
    box_bg ='#A5D6A7', 
    info_box_bg = '#D8DEE9'
  )
)
##########################


header <- dashboardHeader(title = "Demand Forecast")

sidebar <- dashboardSidebar(
    
    # load data
    fileInput('file1', label=h6('Upload Data File'),
              accept=c('text/csv', 
                       'text/comma-separated-values,text/plain', 
                       '.csv')),

    tags$hr(),
    
    radioButtons("Diffusionmodel", label = h6("Choose the Diffusion model"),
                 choices = list("Bass" = 'bass', "Gumbel" = 'gumbel',
                                "Logistic" = 'logis', "Exponential" = 'exp'),selected = 'bass'),
    
    radioButtons("Estimationmethod", label = h6("Choose the Estimation method"),
                 choices = list("OLS" = 'ols', "Q-Q plot" = 'qqplot',
                                "MLE" = 'mle'), selected = 'ols'),

    checkboxInput("weekend_effect",
              label = h6("Weekend adjustment"),
              value = FALSE),

    
    numericInput("bins",
                 label = h6("Forecast Range:"),
                 value = 5,  # 기본값
                 min = 1,    # 최소값
                 max = 30    # 최대값
    ),
    p("Enter a value between 1 and 30")  # 설명 추가
    
)

body <- dashboardBody(
  use_theme(mytheme),
    fluidRow(
        box(title = 'Time-Series plot', status = 'primary', width = 6,
            plotOutput("Plot1")),
        box(title = 'Cumulative Time-Series Plot', status = 'primary', width = 6,
            plotOutput("Plot2"))
        ),
    fluidRow(
        box(title = "Predict values", status = "primary", width = 3,
          tableOutput('predict')),
  valueBoxOutput("m"),
  valueBoxOutput("p"),
  valueBoxOutput("q")
    ))

ui <- dashboardPage(header, sidebar, body)


### estimate parameter
# a,b,c -> m, p, q
bass.coef <- function(coef){
    a <- coef[[1]]; b <- coef[[2]]; c <- coef[[3]]
    m <- (-b-sqrt(b^2-4*a*c))/(2*c)
    p <- a/m
    q <- b+p
    return(list(p=p, q=q, m=m))
}

logit.coef <- function(coef){
    a <- coef[[1]]; b <- coef[[2]]
    q <- a
    m <- -a/b
    return(list(p=0, q=q, m=m))
}

gumbel.coef <- function(coef){
    a <- coef[[1]]; b <- coef[[2]]
    q <- -b
    m <- exp(-a/b)
    return(list(p=0, q=q, m=m))
}

exp.coef <- function(coef){
    a <- coef[[1]]; b <- coef[[2]]
    p <- -b
    m <- -a/b
    return(list(p=p, q=0, m=m))
}


bass.qq.coef <- function(coef){
    m <- coef[[1]]
    k<- coef[[2]]; c<- coef[[3]]
    p <- k/(1+c)
    q <- c*p
    return(list(m=m, p=p, q=q,mu=mu,sigma=sigma))  ##mu=k, sig=c
}

logit.qq.coef <- function(coef){
    m <- coef[[1]]
    mu <- coef[[2]]; sigma <- coef[[3]]
    q <- 1/sigma
    return(list(m=m, p=0, q=q,mu=mu,sigma=sigma))
}

gumbel.qq.coef <- function(coef){
    m <- coef[[1]]
    mu <- coef[[2]]; sigma <- coef[[3]]
    q <- 1/sigma
    return(list(m=m, p=0, q=q,mu=mu,sigma=sigma))
}


exp.qq.coef <- function(coef){
    m <- coef[[1]]
    mu <- coef[[2]]; sigma <- coef[[3]]
    p<- 1/sigma
    return(list(m=m, p=p, q=0,mu=mu,sigma=sigma))
}



# OLS
ols <- function(data, model){
    if(model=="bass"){
        lm.fit <- lm(St ~ y_lag + I(y_lag^2), data)
        par <- bass.coef(lm.fit$coef)
    }
    else if(model=="logis"){
        lm.fit <- lm(St ~ y_lag + I(y_lag^2) -1, data)
        par <- logit.coef(lm.fit$coef)
    }
    else if(model=="gumbel"){
        data <- data[data$y_lag>0,]
        lm.fit <- lm(St ~ y_lag + I(y_lag*log(y_lag))-1, data)
        par <- gumbel.coef(lm.fit$coef)
    }
    else if(model=="exp"){
        lm.fit <- lm(St ~ y_lag, data)
        par <- exp.coef(lm.fit$coef)
    }
    return(par)
    
}



# q-q plot
findm_gum <- function(x,r,m){
    u <- r/(m+1)
    g <- -log(-log(u))
    fit <- lm(x~g)
    return(summary(fit)$r.squared)
}

findm_logis <- function(x,r,m){
    u <- r/(m+1)
    g <- log(u/(1-u))
    fit <- lm(x~g)
    return(summary(fit)$r.squared)
}

findm_exp <- function(x,r,m){
    u <- r/(m+1)
    g <- -log((1-u))
    fit <- lm(x~g)
    return(summary(fit)$r.squared)
}


qqplot <- function(data, model){
    if (model=="bass"){
        lm.fit <- lm(St ~ y_lag + I(y_lag^2), data)
        a<-summary(lm.fit)$coef[1]
        b<-summary(lm.fit)$coef[2]
        c <- summary(lm.fit)$coef[3]
        if(c>=0){
            return(NA)}
        m <- (-b-sqrt(b^2-4*a*c))/(2*c)
        p <- a/m
        q <- -c*m
        
        #data$Ur<-data$Yt/(m+1)
        qq.fit <- nls(t ~ log((1+c*(Yt/(m+1)))/(1-(Yt/(m+1))))/k,
                      data,start=list(m=m,k=p+q,c=q/p))
        par<-bass.qq.coef(c(summary(qq.fit)$coef))
    }
    
    
    else if(model=="logis"){
        lm.fit <- lm(St ~ y_lag + I(y_lag^2) -1, data)
        a<-summary(lm.fit)$coef[1]
        b<-summary(lm.fit)$coef[2]
        m.first<- -a/b
        
        r <- c()
        seq<-seq(m.first-10000,m.first+10000,1000)
        
        for(i in seq){
            r <- c(r,findm_logis(data$t, data$Yt, i))
        }
        
        m <- seq[which.max(r)]
        
        Ur<-data$Yt/(m+1)
        data$logis<-log(Ur/(1-Ur))
        qq.fit <- lm(t~logis, data)
        par<-logit.qq.coef(c(m,qq.fit$coef))
    }
    else if(model=="gumbel"){
        data <- data[data$y_lag>0,]
        lm.fit <- lm(St ~ y_lag + I(y_lag*log(y_lag))-1, data)
        a<-summary(lm.fit)$coef[1]
        b<-summary(lm.fit)$coef[2]
        m.first<-exp(a/(-b))
        
        r <- c()
        seq<-seq(m.first-10000,m.first+10000,1000)
        
        for(i in seq){
            r <- c(r,findm_gum(data$t, data$Yt, i))
        }
        
        m <- seq[which.max(r)]
        #q <- -b
        
        
        Ur<-data$Yt/(m+1)
        data$gumbel<--log(-log(Ur))
        qq.fit <- lm(t ~ gumbel, data)
        par <- gumbel.qq.coef(c(m,qq.fit$coef))
    }
    else if(model=="exp"){
        lm.fit <- lm(St ~ y_lag, data)
        a<-summary(lm.fit)$coef[1]
        b<-summary(lm.fit)$coef[2]
        m.first<-a/(-b)
        
        r <- c()
        seq<-seq(m.first-10000,m.first+10000,1000)
        
        for(i in seq){
            r <- c(r,findm_exp(data$t, data$Yt, i))
        }
        m <- seq[which.max(r)]
        
        Ur<-data$Yt/(m+1)
        data$exp<--log(1-Ur)
        qq.fit <- lm(t ~ exp, data)
        par <- exp.qq.coef(c(m,qq.fit$coef))
    }
    return(par)
    
}



# mle
mle <-function(data, model){
    initial <- c(qqplot(data, model)[c("m","mu","sigma")])
    
    loglik <-function(par){
        m<-par[1]; mu<-par[2]; sigma<-par[3]
        max <- nrow(data)
        if(model=='logis'){
            data["Ft"] <- exp((data$t-mu)/sigma)/(1+exp((data$t-mu)/sigma))
            data[1,'prob'] <- data[1,'Ft'] - exp((0-mu)/sigma)/(1+exp((0-mu)/sigma))
        }
        
        else if(model=='gumbel'){
            data["Ft"] <- exp(-exp(-(data$t-mu)/sigma))
            data[1,"prob"] <- data[1,"Ft"]- exp(-exp(-(0-mu)/sigma))
        }
        
        
        data[2:max,'prob'] <- diff(data$Ft)
        data[(max+1),'prob'] <- 1-data[max, 'Ft']
        data[(max+1), 'St'] <- m-sum(data$St[1:max])
        data <- data %>% mutate(v1 = lfactorial(St), v2=St*log(prob))
        loglik <- lfactorial(m) - sum(data$v1) + sum(data$v2)
        return(-loglik)
    }
    
    
    if(model=="logis"){
        mle <- optim(par=initial, loglik, hessian=T, method= "BFGS")
        par <- mle$par
        res = logit.qq.coef(par)
    }
    else if(model=="gumbel"){
        mle <- optim(par=initial, loglik, hessian=T, method= "BFGS")
        par <- mle$par
        res = gumbel.qq.coef(par)
    }
    
    return(res)
}


# result of parameter predict 
result <- function(data, method, model){
    if(method=="ols"){
        res <- ols(data, model) #list(m,p,q)
    }
    else if(method=="qqplot"){
        res <- qqplot(data, model) #list(m,p,q,mu,sigma)
    }
    else if(method=="mle"){
        res <- mle(data, model) #list(m,p,q,mu,sigma)
    }
    return(res)
    
}


# predict
pred <- function(par, model, step, n){ # n=nrow(data)
    mu <- par$mu; sigma <- par$sigma; m <- par$m
    pred_res <- data.frame(step=c(1:step)) %>% mutate(t = step+n)
    

    if(model=="logis"){
        pred_res <- pred_res %>% mutate(g = exp(-(t-mu)/sigma)/(1+exp(-(t-mu)/sigma))^2)
    }
    else if(model=="gumbel"){
        pred_res <- pred_res %>% mutate(g = exp(-(t-mu)/sigma)*exp(-exp(-(t-mu)/sigma)))
    }
    else if(model=="exp"){
        pred_res <- pred_res %>% mutate(g = exp(-(t-mu)/sigma))
    }
    pred_res["predict"] <- m/sigma*pred_res["g"]
    
    return(pred_res[c("step","predict")])
}



## server.R ##
server <- function(input, output, session){
  data <- reactive({
    infile <- input$file1
    
    if(is.null(infile)){
      return(NULL)
    }
    
    df <- read.csv(infile$datapath)
    colnames(df) <- c("date", "St")
    df$St <- as.numeric(gsub(",", "", df$St))
    df["Yt"] <- cumsum(df$St)
    df["y_lag"] <- lag(df$Yt, default = 0)
    df["t"] <- order(df$date)
    
    ############################holiday #################################
    
    # 2024년의 시작과 끝 날짜 설정
    start_date <- as.Date(df$date[1])
    end_date <- as.Date(tail(df$date, 1))
    
    # 모든 날짜 생성
    all_dates <- seq.Date(from = start_date, to = end_date, by = "day")
    
    # 한국 공휴일 목록
    holidays <- c("01-01", "03-01", "05-01", "05-05", "06-06", "08-15", "10-03", "10-09", "12-25") 
    
    # df 생성
    holiday <- data.frame(date = all_dates,
                          is_weekend = (weekdays(all_dates) %in% c("토요일", "일요일")), # 주말 판단
                          holidays = format(all_dates, "%m-%d") %in% holidays) # 공휴일 판단
    
    # 대체 공휴일 확인 (일요일과 겹칠 경우만)
    holiday$substitute_holiday <- ifelse(holiday$holidays & weekdays(all_dates) %in% c("일요일"), TRUE, FALSE)
    holiday1 <- ifelse(lag(holiday$substitute_holiday) & lag(holiday$holidays), TRUE, FALSE)
    holiday1[1] <- FALSE
    holiday$isHoliday <- holiday$is_weekend | holiday$holidays | holiday1
    
    # df 정리
    holiday <- holiday %>% select(date, isHoliday) %>% filter(isHoliday == TRUE)
    
    df$date <- as.Date(df$date)    
    df <- df %>% left_join(holiday, by = "date")
    df$isHoliday[is.na(df$isHoliday)] <- 0
    
    df$wday <- wday(df$date)
    df$isHoliday[df$wday %in% c(1, 7)] <- 1
    
    if (input$weekend_effect == TRUE) {
      df$St <- ifelse(df$isHoliday == 1, df$St * 2/3, df$St)  
    }
    
    df["Yt"] <- cumsum(df$St)
    df["y_lag"] <- lag(df$Yt, default = 0)
    
    return(df)
  })
  
  observeEvent(input$Estimationmethod, {
    
    output$predict <- renderTable({
      if(input$Estimationmethod == "ols") {
        return(NULL)
      }
      
      data.df <- data()
      par <- result(data.df, input$Estimationmethod, input$Diffusionmodel)
      pred.df <- pred(par, input$Diffusionmodel, input$bins, nrow(data.df))
    })
    
    output$Plot1 <- renderPlot({
      data.df <- data()
      gg <- ggplot(data.df) + geom_line(aes(x = t, y = St), lwd = 1, col ="black") +
        xlab(input$x) + ylab(input$y) + ggtitle(paste(min(data.df$date), "~", max(data.df$date)))
      
      print(gg)
    })
    
    output$Plot2 <- renderPlot({
      data.df <- data()
      gg <- ggplot(data.df) + geom_line(aes(x = t, y = Yt), lwd = 1, col = "brown3") +
        xlab(input$x) + ylab(input$y) + ggtitle(paste(min(data.df$date), "~", max(data.df$date)))
      
      print(gg)
    })
    
    output$m <- renderValueBox({
      data.df <- data()
      valueBox(
        prettyNum(round(result(data.df, input$Estimationmethod, input$Diffusionmodel)$m, 0), big.mark = ","), "Market potential", icon = icon("user-friends"),
        color = "red"
      )
    })
    
    output$p <- renderValueBox({
      data.df <- data()
      valueBox(
        round(result(data.df, input$Estimationmethod, input$Diffusionmodel)$p, 3), "Coefficient of Innovation", icon = icon("signal"),
        color = "yellow"
      )
    })
    
    output$q <- renderValueBox({
      data.df <- data()
      valueBox(
        round(result(data.df, input$Estimationmethod, input$Diffusionmodel)$q, 3), "Coefficient of Imitation", icon = icon("signal"),
        color = "yellow"
      )
    })
  })
}


shinyApp(ui=ui, server = server)