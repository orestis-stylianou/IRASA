irasa <- function (signal,sampling_rate,segment_size=0.9,segments_number=15){
    library(signal)
    library(matrixStats)
    library(RSEIS)
    library(pracma)
    library(eegkit)
    signal <- detrend(signal)
    h <- seq(from=1.1, to=1.9, by=0.05)
    Ntotal <- length(signal)
    Ndata <- 2^floor(log2(Ntotal*segment_size))
    L <- floor((Ntotal-Ndata)/(segments_number-1))
    nfft <- 2^nextpow2(ceiling(max(h))*Ndata)
    Nfrac <- nfft/2 + 1
    frequencies <- seq(0,sampling_rate/2,length.out=Nfrac)
    
    Smixd <- rep(0,Nfrac) 
    taper <- hanning(Ndata)
    for (segment in 0:(segments_number-1)){
        start <- L*segment+1
        time_series <- signal[start:(start+Ndata-1)]
        if (nfft >= length(time_series)){
            time_series <- c(time_series,rep(0,nfft-length(time_series)))
        }else {
            time_series <- time_series[1:nfft]
        }
        xdft <- fft(time_series*taper)
        xdft <- xdft/(min(nfft,length(time_series)))
        xdft[2:length(xdft)] <- xdft[2:length(xdft)]*2
        Smixd <- Smixd+abs(xdft[1:Nfrac])^2
    }
    
    Smixd <- Smixd/segments_number
    filtered_signal <- eegfilter(signal, sampling_rate, .Machine$double.xmin, sampling_rate/(2*ceiling(max(h))), method = "butter",
                                 order = 3, forwardreverse = TRUE,scale = FALSE, plot = FALSE)
    Sfrac <- matrix(0,nrow=Nfrac,ncol=length(h))
    for (i in 1:length(h)){
        factor <- h[i]
        S_up <- rep(0,Nfrac)
        S_down <- rep(0,Nfrac)
        for (segment in 0:(segments_number-1)){
            start <- L*segment+1
            original <- signal[start:(start+Ndata-1)]
            upsampled <- spline(original,n=length(original)*factor)$y
            if (nfft >= length(upsampled)){
                upsampled <- c(upsampled,rep(0,nfft-length(upsampled)))
            }else {
                upsampled <- upsampled[1:nfft]
            }
            taper <- hanning(length(upsampled))
            xdft_up <- fft(upsampled*taper)
            xdft_up <- xdft_up/(min(nfft,length(upsampled)))
            xdft_up[2:length(xdft_up)] <- xdft_up[2:length(xdft_up)]*2
            S_up <- S_up+abs(xdft_up[1:Nfrac])^2
            
            
            original <- filtered_signal[start:(start+Ndata-1)]
            downsampled <- spline(original,n=length(original)*(1/factor))$y
            if (nfft >= length(downsampled)){
                downsampled <- c(downsampled,rep(0,nfft-length(downsampled)))
            }else {
                downsampled <- downsampled[1:nfft]
            }
            taper <- hanning(length(downsampled))
            xdft_down <- fft(downsampled*taper)
            xdft_down <- xdft_down/(min(nfft,length(downsampled)))
            xdft_down[2:length(xdft_down)] <- xdft_down[2:length(xdft_down)]*2
            S_down <- S_down+abs(xdft_down[1:Nfrac])^2
            
            
        }
        
        S_up <- S_up/segments_number
        S_down <- S_down/segments_number
        Sfrac[,i] <- sqrt(S_up*S_down)
    }
    frequencies <- seq(0,sampling_rate/2,length.out=Nfrac)
    Sfrac <- rowMedians(Sfrac)
    Sosc <- abs(Smixd-Sfrac)
    
    
    return(list(original=Smixd, scale_free=Sfrac, oscillatory=Sosc,frequencies=frequencies))
}
    