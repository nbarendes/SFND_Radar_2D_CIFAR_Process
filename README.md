# Radar Target Generation and Detection


# Goal

In this project:


*   Configure the FMCW waveform based on the system requirements.
*   Define the range and velocity of target and simulate its displacement.
*   For the same simulation loop process the transmit and receive signal to determine the beat signal
*   Perform Range FFT on the received signal to determine the Range
*   Towards the end, perform the CFAR processing on the output of 2nd FFT to display the target.

---

### Frequency Modulated Continous Wave (FMCW)
```
%% FMCW Waveform Generation
%Design the FMCW waveform by giving the specs of each of its parameters.

%Bandwidth (B)
Bsweep = c/2*delta_r;
%Chirp Time (Tchirp) based on the Radar's Max Range
Tchirp = 5.5*(Rmax*2/c);
%Slope (slope) of the FMCW
slope =  Bsweep/Tchirp;


%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells
```
---

### Taget Info:
 ```
%% User Defined Range and Velocity of target

% define the target's initial position and velocity. Note : Velocity
% remains contant
%Initial position
InitialRange = 110;
%Velocity
Velocity = 20;
 ```
 ---

### FFT (1D for the range detection):
* Implement the 1D FFT on the Mixed Signal (Computed from the transimit signal and received signal, from code #line 73 ~ 92)
* Reshape the vector into Nr*Nd array.
* Run the FFT on the beat signal along the range bins dimension (Nr)
* Normalize the FFT output.
* Take the absolute value of that output.
* Keep one half of the signal
* Plot the output
* There should be a peak at the initial position of the target


![Figure 1 FFT output](https://user-images.githubusercontent.com/34095574/88806593-ea164a00-d1b0-11ea-8e0e-63c24ae65111.jpg)



Using the frequency from above, we calculate the range as


```
f = Bsweep*(0:(L/2))/L;
R = (c*Tchirp*f)/(2*Bsweep);
```


Output of Range estimation:


![Figure 2 Range from First FFT](https://user-images.githubusercontent.com/34095574/88806773-221d8d00-d1b1-11ea-9468-5e167e7b6800.jpg)

---

### 2D FTT

The 2D FTT is run on the mixed signal (beat signal) output to generate the Range Doppler Map.

```
sig_fft2 = fft2(Mix,Nr,Nd);
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
```

![Figure 3 surface plot of FFT2](https://user-images.githubusercontent.com/34095574/88806965-66109200-d1b1-11ea-93fe-cdab7c99426b.jpg)

---

### 2D CFAR process

![CA-CFAR](https://user-images.githubusercontent.com/34095574/88807427-ff3fa880-d1b1-11ea-8c32-005e6af3f3d9.jpg)


The main outcome of the coherent processing procedure applied to the received echo signal is the two-dimensional range-Doppler-matrix (RDM), which is the basis for an adaptive constant false alarm rate (CFAR) target detection technique. The well-known one-dimensional cell averaging (CA) CFAR procedure suffers from masking effects in all multitarget situations and is additionally very limited in the proper choice of the reference window length. In contrast the ordered statistic (OS) CFAR is very robust in multitarget situations but requires a high computation power. Therefore two-dimensional CFAR procedure based on a combination of OS and CA CFAR is proposed.


![doppler_range](https://user-images.githubusercontent.com/34095574/88807663-4b8ae880-d1b2-11ea-8f0c-02841e744b8a.jpg)

 The false alarm issue can be resolved by implementing the constant false alarm rate. CFAR varies the detection threshold based on the vehicle surroundings. The CFAR technique estimates the level of interference in radar range and doppler cells “Training Cells” on either or both the side of the “Cell Under Test”. The estimate is then used to decide if the target is in the Cell Under Test (CUT).

The process loops across all the range cells and decides the presence of target based on the noise estimate.The basis of the process is that when noise is present, the cells around the cell of interest will contain a good estimate of the noise, i.e. it assumes that the noise or interference is spatially or temporarily homogeneous. Theoretically it will produce a constant false alarm rate, which is independent of the noise or clutter level



#### Implementation Steps



*   The number of Training Cells in both the dimensions are selected

    **Tr = 10**


    **Td = 8**

*   The number of Guard Cells in both dimensions around the Cell under test (CUT) for accurate estimation are selected

    **Gr = 4**


    **Gd = 4**

*   The offset the threshold by SNR value in dB

    **Offset = 6**

*   Slide the Cell Under Test (CUT) across the complete matrix making sure the CUT has margin for Training and Guard cells from the edges.
*   For each iteration we sum the signal level across the training cells and then average is calculated.
*   Offset is added to the threshold to calculate the new threshold.
*   Next, compare the signal under CUT against this threshold.
*   If the CUT signal level is greater than the Threshold, assign a value of 1, else equate it to zero.
*   Since the cell under test are not located at the edges, due to the training cells occupying the edges, we suppress the edges to zero. Any cell value that is neither 1 nor a 0, assign it a zero.


```
num_cells = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1);
signal_cfar = zeros(Nr/2,Nd);

for i = 1:(Nr/2 - (2*Gr+2*Tr))
    for j = 1:(Nd - (2*Gd+2*Td))
        
        s1 = sum(db2pow(RDM(i:i+2*Tr+2*Gr, j:j+2*Td+2*Gd)),'all');
        s2 = sum(db2pow(RDM(i+Tr:i+Tr+2*Gr, j+Td:j+Td+2*Gd)),'all');    
        noise_level = s1 - s2;
        
        threshold = noise_level/num_cells;      
        threshold = pow2db(threshold) + offset;
        threshold = db2pow(threshold);
        
        signal = db2pow(RDM(i+Tr+Gr, j+Td+Gd));
        
        if (signal <= threshold)
            signal = 0;
        else 
            signal = 1;
        end
        
        signal_cfar(i+Tr+Gr,j+Td+Gd) = signal;  
        
    end
end

```


![Figure 4 CA-CFAR Filtered RDM](https://user-images.githubusercontent.com/34095574/88807782-72e1b580-d1b2-11ea-9b0a-280435ca1e21.jpg)

