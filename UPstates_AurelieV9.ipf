#pragma rtGlobals=1		// Use modern global access method.
Menu "Macros"

	"Analyse All UP states", StimUpStatesNoise3()
	"Extract spike rate from Analysed Up states", ExtractUpParameters()
	"Concatenate Waves", ConcatWaves()
	"Filter Waves",  FilterWaves () 
	"Detect ON_OFF states", detectONandOFFstates()
	"RMS-based UP-DOWN state detection",  detectUPandDOWNstates()
End




Function StimUpStatesNoise3()
	// implements faster method for calcuating running standard deviation
	//additional use of Vm histogram to calculate z scores for noise & Vm
	
	Variable win = 200
	Variable DOWNmin = 500
	Variable UPmin = 800
	Variable NoiseMethods = 1
	string extension = "conc"
	variable UPonsetThresh = 1 
	variable UPOffsetThresh = 0 
	Prompt win,"window size for detection (ms)"
	Prompt DOWNmin, "minimum Down state duration (ms)"
	Prompt UPmin, "minimum Up state duration (ms)"
	Prompt NoiseMethods, "Select method for defining Down state", popup "Current recording; Previous Analysis"
	Prompt extension, "enter string to identify waves to be analysed" 
	Prompt UPonsetThresh, "for UPstate ONSET detection, change threshold" 
	prompt UPOffsetThresh,  "for UPstate OFFSET detection, change threshold (make >0 if Upstates are too long)" 
	DoPrompt "Enter paramaters",  win, DOWNmin, UPmin, NoiseMethods, extension, UPonsetThresh, UPOffsetThresh
	
	DOWNmin/=1000
	Upmin/=1000
	extension = "*" + extension + "*"
	
	//set global strings & variables for wave scrolling
	String/g AnalysisWaves =  WaveList(extension,";","")
	String/g AnalysisWaves1
	
	Variable/g DisplayWaveNo
	Variable/g CurrentWave
	
	//Create display waves
	CurrentWave = 0
	Wave wav = $StringFromList(CurrentWave, AnalysisWaves)
	Duplicate/o wav, DisplayWave
	
	//Create variables for Up state detection
	Variable i = 0
	Variable n
	Variable  down
	Variable up
	Variable stop = 0
	Variable dx = deltax(DisplayWave)
	Variable len = numpnts(DisplayWave)
	Variable point_window = floor(win/1000/dx)
	Variable halfwindow = floor(point_window/2)
	Variable sx1, sx2
	Variable/g start
	Variable/g baseline 
	Variable/g width
	String savedDataFolder
	variable allManual = 0 
	
	
	//Edit display wave list
	DisplayWaveNo = 2	
	if(NoiseMethods==1)
		//Display waves
		Display/w=(220,50, 550, 550) DisplayWave
		ShowInfo
		RenameWindow $S_name, set_cursors
		//Create panel for user definition of DOWN state noise
		NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
		DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
		AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
		DrawText 21,20,"Adjust the cursors to mark DOWN state"
		DrawText 60,55,"Scroll through waves"
		DrawText 40,75," to find representative DOWN state"
		DrawText 71,134,"then press Continue."
		Button button0,pos={84,140},size={92,20},title="Continue"
		Button button0,proc=UserCursorAdjust_ContButtProc1
		Button button2,pos={90,80},size={30,20},title="<<"
		Button button2,proc=UserScrollBack_ContButtonProc1
		Button button3,pos={150,80},size={30,20},title=">>"
		Button button3,proc=UserScrollFwd_ContButtonProc1
		PauseForUser tmp_PauseforCursor,set_cursors
	
		//define UP state criteria using user-set (see above) cursors A and B 
		Duplicate/o DisplayWave, temp, tempsq, temp_sd
		Smooth/M=0 50, temp
		tempsq = temp^2
		Wavestats/q/r=[(pcsr(A)-1-halfwindow), (pcsr(A)-1+halfwindow)] temp
		sx1 = V_avg*V_npnts
		Wavestats/q/r=[(pcsr(A)-1-halfwindow), (pcsr(A)-1+halfwindow)] tempsq
		sx2 = V_avg*V_npnts
		for(n = pcsr(A); n<pcsr(B); n+=1)
			sx1 += temp[n+halfwindow+1]
			sx1 -= temp[n-halfwindow]
			sx2 += tempsq[n+halfwindow+1]
			sx2 -= tempsq[n-halfwindow]
			temp_sd[n] = sqrt( ( (V_npnts * sx2) - (sx1 ^2) )/ (V_npnts * (V_npnts - 1)) )
		endfor	
		Wavestats/q/r=(xcsr(A), xcsr(B)) temp_sd // take average SD from user defined "baseline area"
		baseline = V_avg // take average SD from user defined "baseline timespan"
		width = V_sdev 
		Killwaves/z temp, tempsq, temp_sd
		DoWindow/k set_cursors

		//create folder with detection criteria that can be used for subsequent analysis
		savedDataFolder = GetDataFolder(1)
		setDataFolder root:									
		NewDataFolder /S/O AnalysisSettings
		Variable/g G_Start = start
		Variable/g G_baseline = baseline
		Variable/g G_width = width
		SetDataFolder savedDataFolder	
	else
		savedDataFolder = GetDataFolder(1)
		setDataFolder root:AnalysisSettings
		NVAR G_start, G_baseline,  G_width
		start = G_Start 
		baseline = G_baseline
		width = G_width
		SetDataFolder savedDataFolder
	endif // relates to if NoiseMethods = 1 
        
	prompt allManual, " if you wanna analyse it all by hand, set this to 1" 
	doprompt " if you wanna analyse it all by hand, set this to 1" , allManual		
	Make/o/n=(ItemsInList(AnalysisWaves)) Down_Vm_mean, Down_Vm_width	
	for(i=0; i<ItemsInList(AnalysisWaves); i+=1)
		Wave wav = $StringFromList(i, AnalysisWaves)            
		Duplicate/o wav, $("sd_"+StringFromList(i, AnalysisWaves)), $("Ups_"+StringFromList(i, AnalysisWaves))
		Duplicate/o wav, UP_threshold, UP_noise
		Wave UP_marker = $("Ups_"+StringFromList(i, AnalysisWaves) )
		Wave sd = $("sd_"+StringFromList(i, AnalysisWaves))
		UP_marker = 0
           
		if (allManual ==1)
			UP_marker[50000, 60000] = 1 
		else
           
			//Calculate running noise
			Duplicate/o wav, wavef
			Smooth/M=0 50, wavef
			Duplicate/o wavef, wavefsq
			Wave wavefsq
			wavefsq = wavef^2
			Wavestats/q/r=[(ceil(point_window/2)-halfwindow), (ceil(point_window/2)+halfwindow)] wavef
			sx1 = V_avg*V_npnts
			Wavestats/q/r=[(ceil(point_window/2)-halfwindow), (ceil(point_window/2)+halfwindow)] wavefsq
			sx2 = V_avg*V_npnts		
			for(n = (ceil(point_window/2)+1); n<(len-ceil(point_window/2)); n+=1)
				sx1 += wavef[n+halfwindow+1]
				sx1 -= wavef[n-halfwindow]
				sx2 += wavefsq[n+halfwindow+1]
				sx2 -= wavefsq[n-halfwindow]	
				sd[n] = sqrt( ( (V_npnts * sx2) - (sx1 ^2) )/ (V_npnts * (V_npnts - 1)) )		
			endfor
			sd[0, ceil(point_window/2)] = sd[ceil(point_window/2)+1] 
			sd[(len-floor(point_window/2)),] = sd[(len-floor(point_window/2))-1]
			sd-=baseline // substract baseline SD from current one 
			sd/=width   // divide by the SD of the baseline sd estimate                  
			UP_threshold = (sd>UPonsetThresh)?(1):(0) // this trace will be used to find UPstate ONSET 
			UP_noise = (sd>UPOffsetThresh)?(1):(0) // question mark means if "what's before ==1" then =1, : means "else" 0 
			UP_noise = (sd>UPonsetThresh)?(1):(UP_noise)
		
			//find Up states 

			start=0
			variable counter = 0 
			do 		
				FindLevel/edge=1/q/r=(start,) UP_threshold, 0.8 // if 
				up = V_LevelX
				if(V_flag==1)
					break
				else
					FindLevel/edge=2/q/r=(up,) UP_noise, 0.8
					if(V_flag==1)
						break
					else
						down = V_LevelX
					endif           	
					do
						FindLevel/q/r=(down,down+DOWNmin) UP_threshold, 0.8
						if(V_flag==0)
							FindLevel/q/r=(V_LevelX+0.001,) UP_noise, 0.8
							if(V_flag==1)
								stop=1
							else
								down = V_LevelX
							endif
						else
							stop = 1
						endif
					while(stop==0)
					stop=0
				endif 
				FindLevel/edge=1/q/r=(up, start)  UP_noise, 0.8
				if(V_flag==0)
					up = V_levelX
					if((down-up)<UPmin)
						UP_marker[x2pnt(wav, up), x2pnt(wav, down)] = 0.4
					else
						UP_marker[x2pnt(wav, up), x2pnt(wav, down)] = 1
					endif
				else
					UP_marker[x2pnt(wav, start), x2pnt(wav, down)] = 0.4
				endif
				start = down
				counter +=1 
			while(1)
			Killwaves/z UP_threshold, UP_noise	
		endif  
	endfor
	
	//allow manual check of Up states
	CurrentWave = 0
	AnalysisWaves1 =  WaveList("*UPs*",";","")
	Wave wav = $StringFromList(CurrentWave, AnalysisWaves)
	Duplicate/o wav, DisplayWave
	Wave UP_marker = $StringFromList(CurrentWave, AnalysisWaves1)
	Duplicate/o UP_marker, DisplayWave1	
	Display/w=(220,50, 550, 550) DisplayWave
	ModifyGraph rgb=(0,0,0)
	AppendToGraph/r DisplayWave1
	SetAxis right 0,1
	ShowInfo
	RenameWindow $S_name, set_cursors

	NewPanel/k=2/w=(139,341, 442, 622) as "Pause for cursor"
	DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
	AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
	DrawText 40,95,"Scroll through waves to check accuracy"
	DrawText 81,174,"then press Finsih."
	Button button0,pos={94,190},size={92,20},title="Finish"
	Button button0,proc=UserCursorAdjust_ContButtProc1
	Button button1,pos={8,108},size={282,20},title="Manually define UP state"
	Button button1,proc=UserEdit_ContButtonProc3
	Button button4,pos={8,128},size={282,20},title="Manually define DOWN state"
	Button button4,proc=UserEdit_ContButtonProc4
	Button button2,pos={90,58},size={30,20},title="<<"
	Button button2,proc=UserScrollBack_ContButtonProc1
	Button button3,pos={170,58},size={30,20},title=">>"
	Button button3,proc=UserScrollFwd_ContButtonProc1
	PauseForUser tmp_PauseforCursor,set_cursors
	
	DoWindow/k set_cursors
	
	AnalysisWaves1 =  WaveList("*Ups*",";","")
	//update UP state characteristics measurements
	make/o/n=0   UP_sweep 
	Variable baseline1	
	variable counter47 
	for(n = 0; n<ItemsInList(AnalysisWaves1); n+=1)	
		Wave wav = $StringFromList(n, AnalysisWaves)
		Wave UP_marker = $StringFromList(n, AnalysisWaves1)
		make/o/n=(itemsinlist(AnalysisWaves)) UP_frequency
		Make/O/D/N=0 destWave1
		FindLevels/edge=2/q/d=destwave1  UP_marker, 0.5
		make/o/n=10000 UP_amp, UP_duration, UP_area
		print n 
		duplicate/o wav, filtwav 
		
		Wave W_FindLevels			
		if(V_flag!=2)
			Make/O/D/N=0 coefs;
			FilterIIR/CASC/LO=0.05/ORD=100/COEF coefs, filtwav
			InsertPoints 0,V_LevelsFound, UP_sweep
			UP_sweep[0, (V_LevelsFound-1)] = n
			UP_frequency[n] = V_LevelsFound/60
			counter47  = 0 
			Make/O/D/N=0 destWave
			FindLevels/edge=1/q/d=destWave  UP_marker, 0.5
			do 
				wavestats/q/r=(destWave[counter47], destwave1[counter47]) filtwav
				findlevel/q UP_amp, 0 //find next zero point to insert new value, ugly but whatver
				UP_amp[V_LevelX]  = V_min
				wavestats/q/r=(destWave[counter47]-0.200, destwave[counter47]) filtwav
				UP_amp[V_LevelX]  = V_Avg - UP_amp[V_LevelX]
				UP_duration[V_LevelX] = destwave1[counter47] - destWave[counter47]
				UP_area[V_LevelX]= area( filtwav, destWave[counter47], destWave1[counter47])
				counter47 +=1 
			while (counter47 < V_LevelsFound) 
			
		endif
	endfor
	findlevel UP_amp,0 
	DeletePoints  V_LevelX, numpnts(UP_amp)-V_LevelX, UP_amp, UP_duration, UP_area
	Sort UP_sweep, UP_sweep
	display/k=1 UP_duration
	ModifyGraph mode=4,marker=19,rgb=(0,0,0)
	Label left "UPstate Duration (seconds)"
	Label bottom "UPstate number"
	display/k=1  UP_Area
	ModifyGraph mode=4,marker=19,rgb=(0,0,0)
	Label left "UPstate Area (V^2)"
	Label bottom "UPstate number"
	display/k=1 UP_amp
	ModifyGraph mode=4,marker=19,rgb=(0,0,0)
	Label left "UPstate Amplitude (Volts)"
	Label bottom "UPstate number"
	display/k=1 UP_frequency
	ModifyGraph mode=4,marker=19,rgb=(0,0,0)
	Label left "Upstate Frequency (Hz)"
	Label bottom "sweep number"
	
end
/////////////////

/////////////////
//Common functions
Function UserCursorAdjust_ContButtProc1(ctrlName) : ButtonControl
	String ctrlName
	DoWindow/K tmp_PauseforCursor		// Kill self
End

Function UserScrollBack_ContButtonProc1(ctrlName) : ButtonControl
	String ctrlName
	NVAR CurrentWave
	NVAR DisplayWaveNo
	SVAR AnalysisWaves
	SVAR AnalysisWaves1
	
	CurrentWave-=1
	
	Wave D_Wave = DisplayWave
	Wave New_D_Wave = $(StringFromList(CurrentWave, AnalysisWaves))
	
	
	if(WaveExists(New_D_Wave)!=0)
		D_Wave = New_D_Wave
		if(DisplayWaveNo>1)
			Wave D_Wave1 = DisplayWave1
			Wave New_D_Wave1 = $(StringFromList(CurrentWave, AnalysisWaves1))
			D_Wave1 =  New_D_Wave1
		endif
	else
		CurrentWave+=1
	endif
	if(DisplayWaveNo==50)
	
	wavestats/Q D_Wave
	D_Wave*= D_Wave
	wavestats/Q D_Wave
	SetAxis  left V_min ,V_max
	endif

end

Function UserScrollFwd_ContButtonProc1(ctrlName) : ButtonControl
	String ctrlName
	NVAR CurrentWave
	NVAR DisplayWaveNo = DisplayWaveNo 
	SVAR AnalysisWaves
	SVAR AnalysisWaves1
	
	CurrentWave+=1
	
	Wave D_Wave = DisplayWave
	Wave New_D_Wave = $(StringFromList(CurrentWave, AnalysisWaves))
	if(WaveExists(New_D_Wave)!=0)
		D_Wave = New_D_Wave
		if(DisplayWaveNo>1)
			
			Wave D_Wave1 = DisplayWave1
			Wave New_D_Wave1 = $(StringFromList(CurrentWave, AnalysisWaves1))
			D_Wave1 =  New_D_Wave1
		endif
	else
		CurrentWave-=1
	endif
	if(DisplayWaveNo==50)
	
	wavestats/Q D_Wave
	D_Wave*=D_Wave
	wavestats/Q D_Wave
	SetAxis  left V_min ,V_max
	endif
end

Function UserEdit_ContButtonProc2(ctrlName) : ButtonControl
	String ctrlName
	NVAR CurrentWave
	NVAR start
	SVAR AnalysisWaves
	SVAR AnalysisWaves1
	
	Wave UPwav = $(StringFromList(CurrentWave, AnalysisWaves))
	Wavestats/q/r=(0, start-0.05) UPwav
	Variable baseline  = V_avg
	
	Wave D_Wave1 = DisplayWave1
	Wave New_D_Wave1 = $(StringFromList(CurrentWave, AnalysisWaves1))
	D_Wave1 = 0; D_Wave1[x2pnt(UPwav, start), x2pnt(UPwav, xcsr(A))] = 1
	New_D_Wave1 = 0; New_D_Wave1[x2pnt(UPwav, start), x2pnt(UPwav, xcsr(A))] = 1
	
	
	Wave UP_duration; UP_duration[CurrentWave] = xcsr(A)-start
	Wave UP_area; UP_area[CurrentWave] = area(UPwav, start, xcsr(A)) - (baseline*(xcsr(A)-start))
	
end

Function UserEdit_ContButtonProc3(ctrlName) : ButtonControl
	String ctrlName
	NVAR CurrentWave
	NVAR start
	SVAR AnalysisWaves
	SVAR AnalysisWaves1
	
	Wave UPwav = $(StringFromList(CurrentWave, AnalysisWaves))
	//Wavestats/q/r=(0, start-0.05) UPwav
	//Variable baseline  = V_avg
	
	Wave D_Wave1 = DisplayWave1
	Wave New_D_Wave1 = $(StringFromList(CurrentWave, AnalysisWaves1))
	D_Wave1[x2pnt(UPwav, xcsr(A)), x2pnt(UPwav, xcsr(B))] = 1
	New_D_Wave1[x2pnt(UPwav, xcsr(A)), x2pnt(UPwav, xcsr(B))] = 1
end

Function UserEdit_ContButtonProc4(ctrlName) : ButtonControl
	String ctrlName
	NVAR CurrentWave
	NVAR start
	SVAR AnalysisWaves
	SVAR AnalysisWaves1
	
	Wave UPwav = $(StringFromList(CurrentWave, AnalysisWaves))
	
	Wave D_Wave1 = DisplayWave1
	Wave New_D_Wave1 = $(StringFromList(CurrentWave, AnalysisWaves1))
	D_Wave1[x2pnt(UPwav, xcsr(A)), x2pnt(UPwav, xcsr(B))] = 0
	New_D_Wave1[x2pnt(UPwav, xcsr(A)), x2pnt(UPwav, xcsr(B))] = 0
	
end


macro ConcatWaves(nrWaves, identifier)
	Variable nrWaves = 30
	string identifier = "filt"
	Prompt nrWaves, "enter number of waves to concatenate"
	Prompt identifier, "enter term to extract waves, i.e. Vm, or filt" 
	
	String wn,wl
	identifier = "*" + identifier + "*"
	variable counter, intraCounter
	String name, concName
	wl = WaveList(identifier, ";", "")
	counter = 0 

	do 
		wn = StringFromList(counter, wl, ";") 
		if (strlen(wn) == 0) 
			break 
		endif
		
		name = "conc_" + num2str(counter/nrWaves) 
		print name
		wn = StringFromList(counter, wl, ";") 
		duplicate/o $wn, $name
		intraCounter = 0
		concName =   wn  	
		do 
			
			intraCounter+=1
			wn = StringFromList(counter + intraCounter, wl, ";") 
			if (strlen(wn) == 0) 
				print " ATTENTION: last group shorter than others" 
				break 
			endif
			concName = concName + ";" + wn 	
		while (intracounter< (nrWaves - 1))
		concName += ";"
		print concName
		concatenate/o/np  concName, $name
	
		name = "conc_" + num2str(counter) 
	
		counter +=nrWaves 
	while (1)
	
	// move stuff to new folder 
	variable Counter2 = 0 
	killdatafolder/Z conc
	newdatafolder conc
	wl = wavelist ("*conc*",";","")
	Counter2 = 0 
	do 
		wn=stringfromlist(Counter2,wl)
		MoveWave $(wn),:conc:
		Counter2 +=1 
	while(Counter2 < itemsinlist(wl))
			
end


Macro FilterWaves (hiPass, loPass, medianSmooth, notch,move) 
	Variable hiPass = 0
	Variable loPass = 0
	Variable medianSmooth = 0
	variable notch =1
	variable move = 0 
	prompt hiPass, "if nonzero sets highpass cutoff frequency in Hz for highpassfilter (e.g. 5)  "
	prompt loPass, "if nonzero sets lowpass cutoff frequency in Hz for lowpass filter (e.g. 400) "
	prompt mediansmooth, "if nonzero, determine number of POINTS (0.1ms) used for smoothing "
	prompt notch, "if nonzero, applies 50 hz notch filter" 
	prompt move, "set to 1 to move filtered data to new folder" 
	
	string extension = "conc" 
	extension = "*" + extension + "*"
	string wn,wl, name 
	variable counter = 0, n 
	hiPass /= 10000
	loPass/=10000

	wl = WaveList(extension, ";", "")
	do 
		wn = stringfromlist(counter,wl)
		if (strlen(wn) == 0) 
			break 
		endif
			
		name = wn + "filt" 
		print name 
		n = numpnts($name)
		
		duplicate/o $wn, $name, helper
		
		if (notch ==0)  // just regression for noise 
			LineNoise($name)
			$name = signalf1 
		endif 
		// phase conserving notch filter 
		if (notch ==1)
			Make/O/D/N=0 coefs; DelayUpdate
			FilterIIR/N={0.005,25}/COEF coefs, $name
			duplicate/o $name, helper2
			helper2 = $name[numpnts(helper2)-1-p]
			FilterIIR/N={0.005,25}/COEF coefs, helper2
			$name = helper2 [numpnts(helper2)-1-p]
		endif 		
		if (hiPass + loPass > 0)
			duplicate/o $name, helper 
			helper = $name[numpnts(helper)-1-p]
			Concatenate/np/o   {$name,helper}, helper2
			duplicate/o helper2, helper3 
			Concatenate/np/o   {helper2,helper3}, helper 
			
		endif
		
		if (hiPass !=0)
			Make/O/D/N=0 coefs; DelayUpdate
			FilterIIR/CASC/HI=(hiPass)/ORD=1/COEF coefs, helper
			duplicate/o helper, helper2 
			n=numpnts(helper)
			helper2 = helper[n-1-p]
			FilterIIR/CASC/HI=(hiPass)/ORD=1/COEF coefs, helper2
			helper = helper2[n-1-p]	
			duplicate/o/r=[(numpnts($wn)*2),(numpnts($wn)*3)] helper, $name 
			setscale/p x,0,deltax($wn), $name
		endif	
	
		if (loPass !=0)
			Make/O/D/N=0 coefs; DelayUpdate
			duplicate/o  helper, helper2 
			FilterIIR/CASC/LO=(loPass)/ORD=20/COEF coefs, helper
			helper2 = helper[numpnts(helper2)-1-p]
			FilterIIR/CASC/LO=(loPass)/ORD=20/COEF coefs, helper2			
			helper = helper2[numpnts(helper2)-1-p]
			duplicate/o/r=[(numpnts($wn)*2),(numpnts($wn)*3)] helper, $name 
			setscale/p x,0,deltax($wn), $name
		endif 				
		if (medianSmooth != 0)
			Smooth/M=0, medianSmooth, $name
		endif 
		counter += 1
	while(counter < itemsinlist(wl))
	
	if (move == 1)
		// move stuff to new folder 
		variable Counter2 = 0 
		killdatafolder/Z filtered
		newdatafolder filtered
		wl = wavelist ("*filt*",";","")
		Counter2 = 0 
		do 
			wn=stringfromlist(Counter2,wl)
			MoveWave $(wn),:filtered:
			Counter2 +=1 
		while(Counter2 < itemsinlist(wl))
	endif
	
endmacro
	


	

#pragma rtGlobals=3		// Use modern global access method and strict wave access.
Function LineNoise(signal)

	Wave signal 
	Duplicate/o signal, sinWave, cosWave, noise, signalf

	Variable i, sinCoef, cosCoef, n
	for(i=49.4; i<=50.6; i+=0.2)
		sinWave = sin(2*pi*x*i)
		cosWave = cos(2*pi*x*i)
	
		StatsLinearRegression/PAIR sinWave, signal
		Wave W_StatsLinearRegression
		sinCoef = W_StatsLinearRegression[2]
	
		StatsLinearRegression/PAIR cosWave, signal
		Wave W_StatsLinearRegression
		cosCoef = W_StatsLinearRegression[2]
	
		noise = sinCoef*sin(2*pi*i*x) + cosCoef*cos(2*pi*x*i)
	
		signalf-=noise
	
	endfor


	Duplicate/o signal, signalf1
	FFT/OUT=1/DEST=signalFFT signal
	signalFFT = r2polar(signalFFT)

	FFT/OUT=5/DEST=tempPhase signal

	FFT/OUT=3/DEST=tempAmp signal
	tempAmp[x2pnt(tempAmp,49.2),  x2pnt(tempAmp,51)] = tempAmp[x2pnt(tempAmp,49.2)] + (x-49.2)*(tempAmp[x2pnt(tempAmp,51)] - tempAmp[x2pnt(tempAmp,49.2)])

	signalFFT = cmplx(tempAmp[p], tempPhase[p])
	signalFFT = p2rect(signalFFT)

	ifft/DEST= signalf1 signalFFT

	Killwaves/z noise, cosWave, sinWave, signalFFT, tempPhase, tempAMp
end
	
	

	
/////////////////

Function ExtractUpParameters()

	//set global strings
	String/g AnalysisWaves =  WaveList("conc*",";","")
	String/g UPWaves =  WaveList("Ups*",";","")

	if (itemsinlist (UPWaves) - itemsinlist (AnalysisWaves ) != 0)
		print "ERROR: uneven number of waves in  wavelists, better go fix it you muppet!!! "
		print "difference=" + num2str(itemsinlist (UPWaves) - itemsinlist (AnalysisWaves ))
		return (0)
	endif 

	//extract spike rate
	variable sweeps = itemsinlist(AnalysisWaves)
	wave UP_duration = UP_duration
	duplicate/o UP_duration, UP_spikerate,UP_spikeCount
	Variable i, start, stop, helper, UpCounter
	helper = 0
	UpCounter = 0 

	for(i=0; i<sweeps; i+=1)
		Wave wav1 = $StringFromList(i, AnalysisWaves)
		Wave UPwav = $StringFromList(i, UPWaves)
		FindLevel/edge=1/q/r=(0,) UPwav, 0.5
		if(V_flag==0)
	
			Duplicate/O wav1,wav; DelayUpdate
			Make/O/D/N=0 coefs; DelayUpdate
			FilterIIR/CASC/LO=0.3/HI=0.03/ORD=100/COEF coefs, wav
			duplicate wav, $("HpFilt_" + num2str(i))

	
			// when first trace with UPstates appears, the user can define DOWNstate noise	
			if (helper == 0)
				helper = 1 
				Duplicate/o wav, DisplayWave
				//Display waves
				Display/w=(220,50, 550, 550) DisplayWave
				appendtograph UPwav
				wavestats/Q DisplayWave
				SetAxis  left (V_avg- 9*V_sdev) , (V_avg+ 9*V_sdev)
				ModifyGraph rgb(DisplayWave)=(0,0,0)
				ShowInfo
				RenameWindow $S_name, set_cursors
				//Create panel for user definition of DOWN state noise
				NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
				DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
				AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
				DrawText 21,20,"Adjust the cursors to mark noise"
				DrawText 60,55,"Scroll through waves"
				DrawText 40,75," to find representative baseline"
				DrawText 71,134,"then press Continue."
				Button button0,pos={84,140},size={92,20},title="Continue"
				Button button0,proc=UserCursorAdjust_ContButtProc1
				Button button2,pos={90,80},size={30,20},title="<<"
				Button button2,proc=UserScrollBack_ContButtonProc1
				Button button3,pos={150,80},size={30,20},title=">>"
				Button button3,proc=UserScrollFwd_ContButtonProc1
				PauseForUser tmp_PauseforCursor,set_cursors
				//define spike threshold
				Wavestats/q/r=(xcsr(A), xcsr(B)) DisplayWave
				Variable/g spike_threshold = -5*V_sdev
				DoWindow/k set_cursors
			endif 
			// end of downstate noise defintion 
			// the below happens to every wave (finding UPstate spike rate) 
			start =V_LevelX // find beginning of upstate 
			FindLevel/edge=2/q/r=(start+0.05,) UPwav, 0.5 // find end of UPstate
			stop = V_levelX
			FindLevels/Q/edge=2/r=( (start + 0.07), stop ) wav, spike_threshold //find all spikes
			UP_spikerate[UpCounter] = V_LevelsFound/(stop-start)
			UP_spikeCount[UpCounter] = V_LevelsFound
			UpCounter +=1 
			//check for second Up state
			do 
				FindLevel/edge=1/q/r=(stop+0.05,) UPwav, 0.5
				if (V_flag==1)
					break
				elseif(V_flag==0) //repeat as above 
					start =V_LevelX
					FindLevel/edge=2/q/r=(start+0.05,) UPwav, 0.5
					stop = V_levelX
					FindLevels/Q/edge=2/r=( (start + 0.07), stop ) wav, spike_threshold //find all spikes
					UP_spikerate[UpCounter] = V_LevelsFound/(stop-start)
					UP_spikeCount[UpCounter] = V_LevelsFound
					UpCounter +=1 
				endif
			while (1)
		endif	
	endfor
end

Function detectONandOFFstates()

	variable UPmin = 300	
	variable DOWNmin = 300
	string extension = "conc" 
	variable displayStuff =1 
	Variable percentile = 0.01
	Variable minSpikeNum = 5 
	Variable MedianOrMode = 0 

	prompt UPmin, "minimum UPduration  in ms"
	prompt DOWNmin, "minimum UPduration  in ms "
	prompt extension, "part of name to identify things "
	prompt displayStuff, "set to 1 to look through results in the end"
	prompt percentile, " Cutoff percentile UPstate detection (i.e. 0.5 = if spikerate < 95% perecentile of poisson then DOWNstate)"
	prompt minSpikeNum, " enter minimum spike rate of an UPstate"
	prompt MedianOrMode, "to exlude false UPstates, use median (0) or Mode (1) of firing rate distribution "

	doPrompt "enter min up and DOWN state duration", UPmin, DOWNmin, extension, displayStuff, percentile, minSpikeNum, MedianOrMode

	string extension2 =  "*"+ extension +  "*"
	String deletewaves =  WaveList(("!"+extension2),",","") 
	print deletewaves 
	killwaves/Z deletewaves

	

	//set global strings panda
	String/g AnalysisWaves =  WaveList(extension2,";","") 
	variable/g CurrentWave
	variable/g DisplayWaveNo = 0 
	string/g AnalysisWaves1
	//create local variables 
	Variable i, start, stop, helper, Counter, yesno, msCounter, potential, UPSpikeRateMean, UPspikeSD, n 
	Variable threshExist = 0 

	helper = 0
	CurrentWave = 0
	String cmd, content 
	AnalysisWaves =  WaveList("*filt",";","") 
	
	// if not yet done, bandpassfilter signal 300-3000 Hz 
	Counter = 0 
	Variable reSampleRateSpikes = 4000
	UPmin*= (reSampleRateSpikes/1000)
	DOWNmin*=(reSampleRateSpikes/1000)
	if (itemsinlist(AnalysisWaves)==0) 	
		AnalysisWaves =  WaveList(extension2,";","") 
		do 
			if (strlen(stringfromlist(Counter, AnalysisWaves)) ==0)
				break
			endif
			duplicate/o $stringfromlist(Counter, AnalysisWaves),  $(stringfromlist(Counter, AnalysisWaves) + "filt")
			resample/rate=(reSampleRateSpikes) $(stringfromlist(Counter, AnalysisWaves) + "filt")
			wave filtWave = $(stringfromlist(Counter, AnalysisWaves) + "filt")
			duplicate/o filtwave, helpWave
			helpWave =filtWave[numpnts(filtWave)-1-p]
			Concatenate/np/o   {filtWave,helpWave}, helpWave2
			duplicate/o helpWave2, helpWave3 
			Concatenate/np/o   {helpWave2,helpWave3}, helpWave 
			Make/O/D/N=0 coefs; DelayUpdate
			FilterIIR/CASC/HI=(300/reSampleRateSpikes)/ORD=12/COEF coefs, helpWave 
			duplicate/o helpWave, helpWave2 
			n=numpnts(helpWave2)
			helpWave2 = helpWave[n-1-p]
			FilterIIR/CASC/HI=(300/reSampleRateSpikes)/ORD=12/COEF coefs, helpWave 
			helpWave = helpWave2[n-1-p]	
			duplicate/o/r=[(numpnts(filtwave)*2),(numpnts(filtwave)*3)] helpWave, filtWave
			setscale/p x,0,(1/reSampleRateSpikes), filtwave

			Counter +=1 	
		while (1) 
		//sprintf cmd, "FilterWaves (300,3000, 0, 0,0)"
		//Execute cmd	
	endif 
	
	// set up list of waves to be analysed
	AnalysisWaves =  WaveList("*filt",";","") 
	variable sweeps = itemsinlist(AnalysisWaves)
	make/o/n=0 OverAllUpRates, OverAllRates
	
	// MAIN ROUTINE 
	do 
		for(i=0; i<sweeps; i+=1)
			print ("currently processing sweep_nr:" + num2str(i))
			AnalysisWaves =  WaveList("*filt",";","") 
			Wave wav1 = $StringFromList(i, AnalysisWaves)	
			Duplicate/O wav1,wav	
			// when first trace , the user needs to define DOWNstate noise	
			if (exists("spike_threshold")!=0 && i ==0 )
				prompt threshExist, "enter 1 to use previously defined spike threshold, else choose 0" 
				doPrompt  "enter 1 to use previously defined spike threshold, else choose 0" , threshExist
			endif
		 
			if (i== 0&&threshExist!=1)
				do 
					Duplicate/o wav, DisplayWave	
					//Display waves
					Display/w=(220,50, 550, 550) DisplayWave
					wavestats/Q DisplayWave
					SetAxis  left (V_avg- 9*V_sdev) , (V_avg+ 9*V_sdev)
					ModifyGraph rgb(DisplayWave)=(0,0,0)
					ShowInfo
					RenameWindow $S_name, set_cursors
					//Create panel for user definition of DOWN state noise
					NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
					DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
					AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
					DrawText 21,20,"Adjust the cursors to mark noise"
					DrawText 60,55,"Scroll through waves"
					DrawText 40,75," to find representative baseline"
					DrawText 71,134,"then press Continue."
					Button button0,pos={84,140},size={92,20},title="Continue"
					Button button0,proc=UserCursorAdjust_ContButtProc1
					Button button2,pos={90,80},size={30,20},title="<<"
					Button button2,proc=UserScrollBack_ContButtonProc1
					Button button3,pos={150,80},size={30,20},title=">>"
					Button button3,proc=UserScrollFwd_ContButtonProc1
					PauseForUser tmp_PauseforCursor,set_cursors
					//define spike threshold
					Wavestats/q/r=(xcsr(A), xcsr(B)) DisplayWave	
					Variable/g spike_threshold = -5*V_sdev
					// allow user to inspect spike threshold 
					DoWindow/k set_cursors
					duplicate/o DisplayWave, thresh 
					thresh = spike_threshold
					CurrentWave = 0
					Display/k=1/w=(220,50, 550, 550)  DisplayWave
					appendtograph thresh
					ModifyGraph rgb(thresh)=(0,0,0)
					RenameWindow $S_name, set_cursors
					NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
					DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
					AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
					DrawText 71,134,"inspect threshold and press continue ."
					Button button0,pos={84,140},size={92,20},title="Continue"
					Button button0,proc=UserCursorAdjust_ContButtProc1
					Button button2,pos={90,80},size={30,20},title="<<"
					Button button2,proc=UserScrollBack_ContButtonProc1
					Button button3,pos={150,80},size={30,20},title=">>"
					Button button3,proc=UserScrollFwd_ContButtonProc1
					PauseForUser tmp_PauseforCursor,set_cursors					
					prompt yesno, "enter 1 to accept threshold" 
					doprompt "enter 1 to accept threshold" , yesno 
					KillWindow set_cursors
					if (yesno == 1)
						break
					endif
				while(1)
			else 
				NVAR spike_threshold =  spike_threshold 
			endif 
		
			// end of downstate noise defintion 
			// the below happens to every wave (finding UPstate spike rate) 	
			duplicate/o wav1, $("Raster_" + num2str(i)), $("UPs_" + num2str(i)), $("UP_rate" + num2str(i)), $("UP_cleaned" + num2str(i))
			wave rasterPlot =  $("Raster_" + num2str(i))		
			wave UPdetect =  $("UPs_" + num2str(i))
			wave UPrate =  $("UP_rate" + num2str(i))	
			wave UPcleaned = $("UP_cleaned" + num2str(i))
			UPdetect = 0 
			UPrate = 0 
			UPcleaned = 0 
		
			// threshold the recording so that all entries above spike threshold are = 1 
			rasterPlot = rasterPlot<spike_threshold 	? 1 : 0
		
			// // Zero all timepoints "during the spike", so that only 1 point per spike remains
			msCounter = 0 
			do 
				FindLevel /EDGE=1/P/Q/R=[msCounter,numpnts(rasterPlot) ] rasterPlot, 0.8 // find next spike 
				if (V_flag == 1) //  if none, stop loop
					break
				endif
				msCounter = V_LevelX // use msCounter to store spike start time 
				FindLevel /EDGE=2/P/Q/R=[V_LevelX,numpnts(rasterPlot) ] rasterPlot, 0.2//find end of spike-like signal
				if (V_flag == 1)// if none, stop loop
					break
				endif
				rasterPlot[msCounter+1, V_LevelX] = 0 // Zero all timepoints "during the spike", so that only 1 point per spike remains
				msCounter = V_LevelX // use msCounter to mark end of last spike, end thus beginning of next search 
			while (1)
		
			// now detect UPstates according to below rule: 
			//1) once a spike is detected and the next DOWNstate(period of no spiking lasting for minDOWN) follows before minUPduration is reached, the whole period is counted as DOWNstate
			msCounter = 0 
			potential = 0 
			do 
				if ((rasterPlot[msCounter] > 0) &&potential == 0) // if there is a spike but we're not in a potential UPstate then mark as potential UPstate
					potential = msCounter 
				elseif (potential !=0) 
					if (sum(rasterPlot,msCounter/reSampleRateSpikes, ((msCounter + DOWNmin)/reSampleRateSpikes))  < 1) // if a downstate occurs 
						if ((msCounter - potential) < UPmin) // if the potential upstate was to short , ditch it 
							potential = 0
						else	 //if it was long enough 
							findlevels/edge=1/q/r=[potential-1,msCounter+1] rasterPlot, 0.8 	// get number of spikes in the UPstate
							if (V_LevelsFound<minSpikeNum) // exclude upstates with too few spikes 
								potential = 0
							print V_LevelsFound
							else
								UPdetect[potential,msCounter] = 1 // if the potential upstate had been long enough: mark it 
								UPrate[potential,msCounter] =  V_LevelsFound/((msCounter - potential)/reSampleRateSpikes) /// save the spike rate (in Hz)
								insertpoints 0,1, OverAllUpRates
								OverAllUpRates[0] = UPrate[potential] 
								potential = 0 
							endif
						endif
					endif 
				endif
				msCounter += 1
			while (msCounter<numpnts(rasterPlot))
			insertpoints 0,1, OverAllRates
			OverAllRates[0] = sum(rasterPlot)/(numpnts(RasterPlot)/reSampleRateSpikes)
		endfor		
		// 2) an UPstate spikerate needs to be at least within the 0.05 percentile of a typical UPstate
		// first allow user to choose a typical upstate 
		AnalysisWaves =  WaveList("*filt",";","") 
		for(i=0; i<sweeps; i+=1)
			wave rasterPlot =  $("Raster_" + num2str(i))		
			wave UPdetect =  $("UPs_" + num2str(i))
			wave UPrate =  $("UP_rate" + num2str(i))	
			wave UPcleaned = $("UP_cleaned" + num2str(i))		
			UPcleaned = UPdetect 
			if (i==0)
				displayStuff = 0 
				prompt displayStuff, " type 1 to manually define 'typlical' upstate (shouldnt be the default choice)" 
				doprompt  "choices...." , displayStuff
				if (displayStuff ==1)		// if user chooses to manually define UPstate
					DisplayWaveNo = 2
					do 	
						duplicate/o UPdetect, DisplayWave1
						AnalysisWaves1 =  wavelist("*UPs*", ";", "") 
						AnalysisWaves=  wavelist("*Raster*", ";", "") 
						wav = rasterPlot 
						Duplicate/o wav, DisplayWave
						//Display waves
						Display/w=(220,50, 550, 550) DisplayWave, DisplayWave1
						wavestats/Q DisplayWave
						SetAxis  left 0, 1
						ModifyGraph rgb(DisplayWave)=(0,0,0)
						ShowInfo
						RenameWindow $S_name, set_cursors
						//Create panel for user definition of UP state 
						NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
						DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
						AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
						DrawText 21,20,"Adjust the cursors to mark a typical UPstate"
						DrawText 60,55,"Scroll through waves"
						DrawText 40,75," to find representative UPstate"
						DrawText 71,134,"then press Continue."
						Button button0,pos={84,140},size={92,20},title="Continue"
						Button button0,proc=UserCursorAdjust_ContButtProc1
						Button button2,pos={90,80},size={30,20},title="<<"
						Button button2,proc=UserScrollBack_ContButtonProc1
						Button button3,pos={150,80},size={30,20},title=">>"
						Button button3,proc=UserScrollFwd_ContButtonProc1
						PauseForUser tmp_PauseforCursor,set_cursors
						//get spike rate 
						findlevels/edge=1/q/r=(xcsr(A), xcsr(B)) wav, 0.8 	// get number of spikes in the UPstate
						if (V_flag != 2) 
							break 
						else 
							print "no spikes found in Upstate, something is no bueno here, buddy! "
						endif 
					while(1) 
					UPSpikeRateMean =  V_LevelsFound/((xcsr(B) - xcsr(A)))
					killwindow set_cursors
				else  // if not user defined UPState
					wavestats/q OverallUpRates 
					if (MedianOrMode == 0)
						UPSpikeRateMean = 	Median(OverAllUpRates, 0, numpnts(OverAllUpRates))
						print "median spike rate is:" + num2str(UPSpikeRateMean) + "_Hz"
					else 
						make/o histo
						Histogram/B=3 OverAllUpRates, histo
						UPSpikeRateMean =  V_maxloc
					endif
				endif
				UPspikeSD = V_avg- 2*V_sdev
				print "standard deviation is:" + num2str(V_sdev) + "_Hz"

				//UPspikeSD=  StatsInvPoissonCDF(  percentile,UPSpikeRateMean )
			endif
			
	
		// now eliminate upstates with spike rate smaller than the entered percentile in the poisson distribution taken from 
			msCounter = 0 
			do 
				FindLevel /EDGE=1/P/Q/R=[msCounter,numpnts(UPdetect)] UPdetect, 0.8 // find next UPstate
				if (V_flag == 1) //  if none, stop loop
					break
				endif
				msCounter = V_LevelX // use msCounter to store spike start time 
				FindLevel /EDGE=2/P/Q/R=[V_LevelX,numpnts(UPdetect) ] UPdetect, 0.2//find end of Upstate
				if (V_flag == 1)// if none, stop loop
					print "no downstate transition found anymore, not sure this should happen" 
					break
				endif
			
				if (UPrate[msCounter+1] < UPSpikeSD) // if upstate spike rate too low, get rid of it 
					UPcleaned[msCounter, V_LevelX] = 0 
					print ("deleted UPstate with " + num2str(UPrate[msCounter+1]) + " Hz spike rate")
				endif 
				msCounter = V_LevelX // use msCounter to mark end of last UPstate,  thus beginning of next search 
			while(1)
			//RasterPlot= RasterPlot < 0.5 ? nan : 1 	
		endfor 
		displayWaveNo = 2
		displayStuff = 0
		prompt displayStuff, " type 1 to examine rasterPlots" 
		doprompt  " type 1 to examine rasterPlots" , displayStuff
		if (displayStuff == 1)
			AnalysisWaves1 = wavelist("*UP_cleaned*", ";", "") 
			AnalysisWaves = wavelist("*raster*", ";", "") 
			Duplicate/o $stringfromlist(0,AnalysisWaves) DisplayWave
			duplicate/o $stringfromlist(0,AnalysisWaves1) DisplayWave1
			//Display waves
			Display/w=(220,50, 550, 550) DisplayWave, DisplayWave1
			wavestats/Q DisplayWave
			SetAxis  left, 0, wavemax(wav1)
			ModifyGraph rgb(DisplayWave)=(43520,43520,43520), rgb(DisplayWave1)=(65280,0,0), lsize(DisplayWave1)=2
			ShowInfo
			RenameWindow $S_name, set_cursors
			//Create panel for user definition of UP state 
			NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
			DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
			AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
			DrawText 60,55,"Scroll through waves"
			DrawText 40,75," to check  algorithm performance"
			DrawText 71,134,"then press Continue."
			Button button0,pos={84,140},size={92,20},title="Continue"
			Button button0,proc=UserCursorAdjust_ContButtProc1
			Button button2,pos={90,80},size={30,20},title="<<"
			Button button2,proc=UserScrollBack_ContButtonProc1
			Button button3,pos={150,80},size={30,20},title=">>"
			Button button3,proc=UserScrollFwd_ContButtonProc1
			PauseForUser tmp_PauseforCursor,set_cursors
			killwindow set_cursors
		endif 	
		
	
		displayStuff = 0
		prompt displayStuff, " type 1 to examine original trace" 
		doprompt  " type 1 to examine original trace" , displayStuff
		Variable count = 0 
		String AnalysisWaves0, wn
		if (displayStuff == 1)
		displayWaveNo = 50 
			AnalysisWaves1 = wavelist("*UP_cleaned*", ";", "") 
			AnalysisWaves0 = wavelist("*conc*", ";", "") 
			AnalysisWaves =  wavelist("this_string_should_not_be_in_any_name", ";", "") 
		     do 
				wn = stringfromlist(count,AnalysisWaves0, ";")
				if (strlen(wn)==0)
					break
				endif 
				if (StringMatch(wn, "*filt*")==0)
					AnalysisWaves  = AddListItem(wn,AnalysisWaves,";",5000)	
				endif
				count+=1 
			while (1)
			Duplicate/o $stringfromlist(0,AnalysisWaves) DisplayWave
			duplicate/o $stringfromlist(0,AnalysisWaves1, ";") DisplayWave1
			//Display waves
			Display/w=(220,50, 550, 550) DisplayWave, DisplayWave1
			wavestats/Q DisplayWave
			Displaywave*= Displaywave
			SetAxis  left, V_min, V_max
			ModifyGraph rgb(DisplayWave)=(43520,43520,43520), rgb(DisplayWave1)=(65280,0,0), lsize(DisplayWave1)=2
			ShowInfo
			RenameWindow $S_name, set_cursors
			//Create panel for user definition of UP state 
			NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
			DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
			AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
			DrawText 60,55,"Scroll through waves"
			DrawText 40,75," to check algorithm performance"
			DrawText 71,134,"then press Continue."
			Button button0,pos={84,140},size={92,20},title="Continue"
			Button button0,proc=UserCursorAdjust_ContButtProc1
			Button button2,pos={90,80},size={30,20},title="<<"
			Button button2,proc=UserScrollBack_ContButtonProc1
			Button button3,pos={150,80},size={30,20},title=">>"
			Button button3,proc=UserScrollFwd_ContButtonProc1
			PauseForUser tmp_PauseforCursor,set_cursors
			killwindow set_cursors
		endif 	
		displayStuff = 0 
		prompt displayStuff, " Happy with result (0) or Run again with different seetings (1)?" 
		doprompt  " inputinput baby" , displayStuff
		if (displayStuff != 1)
			break
		else 
			prompt percentile, " change percentile if you wnat (lower will mean more UPstates will be included)" 
			doprompt  " inputinput baby" , percentile
		endif 
	while(1)
	
	//allow manual check of Up states
	displayStuff = 0 
	prompt displayStuff, " type 1 to manually UPdate UPstate measurements" 
	doprompt  " inputinput baby" , displayStuff
	killwaves/Z filtwav
	if (displayStuff == 1)
		// panda start 
		CurrentWave = 0
		AnalysisWaves1 = wavelist("*UP_cleaned*", ";", "") // this your detection signal (0/1) 
		AnalysisWaves = wavelist("*filt*", ";", "")  // this is your signal
		Wave wav = $StringFromList(CurrentWave, AnalysisWaves)
		Wave UP_marker = $StringFromList(CurrentWave, AnalysisWaves1)
		Duplicate/o $stringfromlist(0,AnalysisWaves) DisplayWave
		duplicate/o $stringfromlist(0,AnalysisWaves1) DisplayWave1	
		Display/w=(220,50, 550, 550) DisplayWave
		ModifyGraph rgb=(0,0,0)
		AppendToGraph/r DisplayWave1
		SetAxis right 0,1
		ShowInfo
		RenameWindow $S_name, set_cursors

		NewPanel/k=2/w=(139,341, 442, 622) as "Pause for cursor"
		DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
		AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
		DrawText 40,95,"Scroll through waves to check accuracy"
		DrawText 81,174,"then press Finsih."
		Button button0,pos={94,190},size={92,20},title="Finish"
		Button button0,proc=UserCursorAdjust_ContButtProc1
		Button button1,pos={8,108},size={282,20},title="Manually define UP state"
		Button button1,proc=UserEdit_ContButtonProc3
		Button button4,pos={8,128},size={282,20},title="Manually define DOWN state"
		Button button4,proc=UserEdit_ContButtonProc4
		Button button2,pos={90,58},size={30,20},title="<<"
		Button button2,proc=UserScrollBack_ContButtonProc1
		Button button3,pos={170,58},size={30,20},title=">>"
		Button button3,proc=UserScrollFwd_ContButtonProc1
		PauseForUser tmp_PauseforCursor,set_cursors
		DoWindow/k set_cursors
	endif 
	
	//update UP state measurements
	// 1) get only original waves into "analysisWaves"
	AnalysisWaves =  wavelist("this_string_should_not_be_in_any_name", ";", "") 
	AnalysisWaves0 = wavelist("*conc*", ";", "") 
	do 
		wn = stringfromlist(count,AnalysisWaves0, ";")
		if (strlen(wn)==0)
			break
		endif 
		if (StringMatch(wn, "*conc*")==0)
			AnalysisWaves  = AddListItem(wn,AnalysisWaves,";",5000)	
		endif
		count+=1 
	while (1)
	 // finished with 1) 
	AnalysisWaves1 = wavelist("*UP_cleaned*", ";", "") 
	
	
	
	
	make/o/n=0   UP_sweep 
	Variable baseline1	
	variable counter47 
	for(n = 0; n<ItemsInList(AnalysisWaves1); n+=1)	
		Wave wav = $StringFromList(n, AnalysisWaves)
		Wave UP_marker = $StringFromList(n, AnalysisWaves1)
		make/o/n=(itemsinlist(AnalysisWaves)) UP_frequency
		Make/O/D/N=0 destWave1
		FindLevels/edge=2/q/d=destwave1  UP_marker, 0.5 // destwave1 will contain downstate transitions
		make/o/n=10000 UP_amp, UP_duration, UP_area
		duplicate/o wav, filtwav 
		
		Wave W_FindLevels			
		if(V_flag!=2)
			Make/O/D/N=0 coefs;
			FilterIIR/CASC/LO=0.05/ORD=100/COEF coefs, filtwav
			InsertPoints 0,V_LevelsFound, UP_sweep
			UP_sweep[0, (V_LevelsFound-1)] = n
			UP_frequency[n] = V_LevelsFound/60
			counter47  = 0 
			Make/O/D/N=0 destWave // // destwave will contain  UPstate downstate transitions
			FindLevels/edge=1/q/d=destWave  UP_marker, 0.5
			do 
				wavestats/q/r=(destWave[counter47], destwave1[counter47]) filtwav
				findlevel/q UP_amp, 0 //find next zero point to insert new value, ugly but whatver
				UP_amp[V_LevelX]  = V_min
				wavestats/q/r=(destWave[counter47]-0.200, destwave[counter47]) filtwav
				UP_amp[V_LevelX]  = V_Avg - UP_amp[V_LevelX]
				UP_duration[V_LevelX] = destwave1[counter47] - destWave[counter47]
				UP_area[V_LevelX]= area( filtwav, destWave[counter47], destWave1[counter47])
				counter47 +=1 
			while (counter47 < V_LevelsFound) 		
		endif
	endfor	
end

Function detectUPandDOWNstates()

	variable UPmin = 300	
	variable DOWNmin = 300
	string extension = "conc" 
	variable displayStuff =1 
	Variable windowSize = 10
	Variable userThresh = 8

	prompt UPmin, "minimum UPduration  in ms"
	prompt DOWNmin, "minimum UPduration  in ms "
	prompt extension, "part of name to identify things "
	prompt displayStuff, "set to 1 to look through results in the end"
	prompt windowSize, "define length of RMS window"
	prompt userThresh, " enter factor (treshold = input*SD)"


	doPrompt "enter min up and DOWN state duration", UPmin, DOWNmin, extension, displayStuff, windowSize, userThresh

	string extension2 =  "*"+ extension +  "*"
	String deletewaves =  WaveList(("!"+extension2),",","") 
	print deletewaves 
	killwaves/Z deletewaves
	variable minSpikeNum = 3 // relatively obsolete from ON OFF detection algorithm 
	

	//set global strings panda
	String/g AnalysisWaves =  WaveList(extension2,";","") 
	variable/g CurrentWave
	variable/g DisplayWaveNo = 0 
	string/g AnalysisWaves1
	//create local variables 
	Variable i, start, stop, helper, Counter, yesno, msCounter, potential, UPSpikeRateMean, UPspikeSD, n, rmsCnt, rmsHelper, m 
	
	Variable threshExist = 0 

	helper = 0
	CurrentWave = 0
	String cmd, content, AnalysisWaves0, wn
	AnalysisWaves =  WaveList("*LOWPASS",";","") 
	
	// if not yet done, lowpassfilter signal 
	Counter = 0 
	Variable newSampleRate = 4000
	UPmin*=(newSampleRate/1000)
	DOWNmin*=(newSampleRate/1000)
	windowSize *=(newSampleRate/1000) // get shit in "pointsclae" relative to new sample rate 
	AnalysisWaves0 = wavelist("*conc*", ";", "") 
	AnalysisWaves =  wavelist("this_string_should_not_be_in_any_name", ";", "") 	
		do 
			wn = stringfromlist(Counter,AnalysisWaves0, ";")
			if (strlen(wn)==0)
				break
			endif 
			if ((StringMatch(wn, "*filt*")==0)&&(StringMatch(wn, "*RMS*")==0) &&(StringMatch(wn, "*LOWPASS*")==0)&&(StringMatch(wn, "*DIF*")==0)&&(StringMatch(wn, "*fit*")==0))
				AnalysisWaves  = AddListItem(wn,AnalysisWaves,";",5000)	
			endif
			Counter+=1 
		while (1)
		Counter = 0
		do // do the filtering 
			if (strlen(stringfromlist(Counter, AnalysisWaves)) ==0)
				break
			endif
			duplicate/o $stringfromlist(Counter, AnalysisWaves),  $(stringfromlist(Counter, AnalysisWaves) + "_LOWPASS")
			resample/rate=(newSampleRate) $(stringfromlist(Counter, AnalysisWaves) + "_LOWPASS")	
			duplicate/o  $(stringfromlist(Counter, AnalysisWaves) + "_LOWPASS"),  $(stringfromlist(Counter, AnalysisWaves) + "_RMS")		
			wave rms = $(stringfromlist(Counter, AnalysisWaves) + "_RMS")
			rms *= rms 
			rmsCnt = 0 
			do // make rolling rms calculatiuons for this wave 
				rmsHelper = sum(rms, pnt2x(rms,rmsCnt), pnt2x(rms,(rmsCnt+windowSize))) 
				rmsHelper/= windowSize
				rmsHelper = sqrt(rmsHelper) 
				rms[rmsCnt, (rmsCnt+windowSize)-1]= rmsHelper
				rmsCnt+=windowSize
				if ((rmsCnt + windowSize) > (numpnts(rms)))
					rmsHelper = sum(rms,pnt2x(rms,rmsCnt), pnt2x(rms,numpnts(rms))) 
					rmsHelper/= windowSize
					rmsHelper = sqrt(rmsHelper) 
					rms[rmsCnt,numpnts(rms)] = rmsHelper
					break
				endif
			while (1)
			MartinFiltFiltIIR ((stringfromlist(Counter, AnalysisWaves) + "_RMS"), newSampleRate, 1,1)
			wave rms_detrended = $(stringfromlist(Counter, AnalysisWaves) + "_RMSfilt")
			rms_detrended = abs(rms_detrended) 
			wavestats/q rms_detrended
			Counter +=1 	
			
		while (1) 
		//sprintf cmd, "FilterWaves (300,3000, 0, 0,0)"
		//Execute cmd	

	// MAIN ROUTINE 
	// set up list of waves to be analysed
	AnalysisWaves =  WaveList("*_RMSfilt",";","") 
	variable sweeps = itemsinlist(AnalysisWaves)	
	do 
		for(i=0; i<sweeps; i+=1)
			print ("currently processing sweep_nr:" + num2str(i))
			Wave wav1 = $StringFromList(i, AnalysisWaves)	
			Duplicate/O wav1,wav	
			// when first trace , the user needs to define DOWNstate noise	
			if (exists("DIF_threshold")!=0 && i ==0 )
				prompt threshExist, "enter 1 to use previously defined DIF threshold, else choose 0" 
				doPrompt  "enter 1 to use previously defined DIF threshold, else choose 0" , threshExist
			endif
		 
			if (i== 0&&threshExist!=1)
				do 
					Duplicate/o wav, DisplayWave	
					//Display waves
					Display/w=(220,50, 550, 550) DisplayWave
					wavestats/Q DisplayWave
					SetAxis  left (V_avg- 9*V_sdev) , (V_avg+ 9*V_sdev)
					ModifyGraph rgb(DisplayWave)=(0,0,0)
					ShowInfo
					RenameWindow $S_name, set_cursors
					//Create panel for user definition of DOWN state noise
					NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
					DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
					AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
					DrawText 21,20,"Adjust the cursors to mark noise"
					DrawText 60,55,"Scroll through waves"
					DrawText 40,75," to find representative baseline"
					DrawText 71,134,"then press Continue."
					Button button0,pos={84,140},size={92,20},title="Continue"
					Button button0,proc=UserCursorAdjust_ContButtProc1
					Button button2,pos={90,80},size={30,20},title="<<"
					Button button2,proc=UserScrollBack_ContButtonProc1
					Button button3,pos={150,80},size={30,20},title=">>"
					Button button3,proc=UserScrollFwd_ContButtonProc1
					PauseForUser tmp_PauseforCursor,set_cursors
					//define spike threshold
					Wavestats/q/r=(xcsr(A), xcsr(B)) DisplayWave	
					Variable/g DIF_threshold = userThresh*V_sdev
					// allow user to inspect spike threshold 
					DoWindow/k set_cursors
					duplicate/o DisplayWave, thresh 
					thresh = DIF_threshold
					CurrentWave = 0
					Display/k=1/w=(220,50, 550, 550)  DisplayWave
					appendtograph thresh
					ModifyGraph rgb(thresh)=(0,0,0)
					RenameWindow $S_name, set_cursors
					NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
					DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
					AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
					DrawText 71,134,"inspect threshold and press continue ."
					Button button0,pos={84,140},size={92,20},title="Continue"
					Button button0,proc=UserCursorAdjust_ContButtProc1
					Button button2,pos={90,80},size={30,20},title="<<"
					Button button2,proc=UserScrollBack_ContButtonProc1
					Button button3,pos={150,80},size={30,20},title=">>"
					Button button3,proc=UserScrollFwd_ContButtonProc1
					PauseForUser tmp_PauseforCursor,set_cursors					
					prompt yesno, "enter 1 to accept threshold" 
					doprompt "enter 1 to accept threshold" , yesno 
					KillWindow set_cursors
					if (yesno == 1)
						break
					endif
				while(1)
			else 
				NVAR DIF_threshold =  DIF_threshold 
			endif 
		
			// end of downstate noise defintion 
			// the below happens to every wave (finding UPstate spike rate) 	
			duplicate/o wav1, $("DIFRaster_" + num2str(i)), $("DIFUPs_" + num2str(i)), $("DIFUP_rate" + num2str(i)), $("DIFUP_cleaned" + num2str(i))
			wave rasterPlot =  $("DIFRaster_" + num2str(i))		
			wave UPdetect =  $("DIFUPs_" + num2str(i))
			wave UPrate =  $("DIFUP_rate" + num2str(i))	
			wave UPcleaned = $("DIFUP_cleaned" + num2str(i))
			UPdetect = 0 
			UPrate = 0 
			UPcleaned = 0 
		
			// threshold the recording so that all entries above spike threshold are = 1 
			rasterPlot = rasterPlot>DIF_threshold	? 1 : 0
		
			// now detect UPstates according to below rule: 
			//1) once a spike is detected and the next DOWNstate(period of no spiking lasting for minDOWN) follows before minUPduration is reached, the whole period is counted as DOWNstate
			msCounter = 0 
			potential = 0 
			do 
				if ((rasterPlot[msCounter] > 0) &&potential == 0) // if there is a spike but we're not in a potential UPstate then mark as potential UPstate
					potential = msCounter 
				elseif (potential !=0) 
					if (sum(rasterPlot,msCounter/newSampleRate, ((msCounter + DOWNmin)/newSampleRate))  < 1) // if a downstate occurs 
						if ((msCounter - potential) < UPmin) // if the potential upstate was to short , ditch it 
							potential = 0
						else	 //if it was long enough 
							findlevels/edge=1/q/r=[potential-1,msCounter+1] rasterPlot, 0.8 	// get number of spikes in the UPstate
							if (V_LevelsFound<minSpikeNum) // exclude upstates with too few spikes 
								potential = 0
							print V_LevelsFound
							else
								UPdetect[potential,msCounter] = 1 // if the potential upstate had been long enough: mark it 
								potential = 0 
							endif
						endif
					endif 
				endif
				msCounter += 1
			while (msCounter<numpnts(rasterPlot))
		endfor		
		
		displayWaveNo = 2
		displayStuff = 0
		prompt displayStuff, " type 1 to examine rasterPlots" 
		doprompt  " type 1 to examine rasterPlots" , displayStuff
		if (displayStuff == 1)
			AnalysisWaves1 = wavelist("*DIFUPs_*", ";", "") 
			AnalysisWaves = wavelist("*raster*", ";", "") 
			Duplicate/o $stringfromlist(0,AnalysisWaves) DisplayWave
			duplicate/o $stringfromlist(0,AnalysisWaves1) DisplayWave1
			//Display waves
			Display/w=(220,50, 550, 550) DisplayWave, DisplayWave1
			wavestats/Q DisplayWave
			SetAxis  left, 0, wavemax(wav1)
			ModifyGraph rgb(DisplayWave)=(43520,43520,43520), rgb(DisplayWave1)=(65280,0,0), lsize(DisplayWave1)=2
			ShowInfo
			RenameWindow $S_name, set_cursors
			//Create panel for user definition of UP state 
			NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
			DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
			AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
			DrawText 60,55,"Scroll through waves"
			DrawText 40,75," to check  algorithm performance"
			DrawText 71,134,"then press Continue."
			Button button0,pos={84,140},size={92,20},title="Continue"
			Button button0,proc=UserCursorAdjust_ContButtProc1
			Button button2,pos={90,80},size={30,20},title="<<"
			Button button2,proc=UserScrollBack_ContButtonProc1
			Button button3,pos={150,80},size={30,20},title=">>"
			Button button3,proc=UserScrollFwd_ContButtonProc1
			PauseForUser tmp_PauseforCursor,set_cursors
			killwindow set_cursors
		endif 	
		
	
		displayStuff = 1
		Variable count = 0 
		if (displayStuff == 1)
		displayWaveNo = 50 
			AnalysisWaves1 = wavelist("*DIFUPs_*", ";", "") 
			AnalysisWaves0 = wavelist("*conc*", ";", "") 
			AnalysisWaves =  wavelist("this_string_should_not_be_in_any_name", ";", "") 
		     do 
				wn = stringfromlist(count,AnalysisWaves0, ";")
				if (strlen(wn)==0)
					break
				endif 
				if (StringMatch(wn, "*filt*")==0&&(StringMatch(wn, "*RMS*")==0) &&(StringMatch(wn, "*LOWPASS*")==0)&&(StringMatch(wn, "*DIF*")==0))
					AnalysisWaves  = AddListItem(wn,AnalysisWaves,";",5000)	
				endif
				count+=1 
			while (1)
			Duplicate/o $stringfromlist(0,AnalysisWaves) DisplayWave
			duplicate/o $stringfromlist(0,AnalysisWaves1, ";") DisplayWave1
			//Display waves
			Display/w=(220,50, 550, 550) DisplayWave, DisplayWave1
			wavestats/Q DisplayWave
			Displaywave*= Displaywave
			SetAxis  left, V_min, V_max
			ModifyGraph rgb(DisplayWave)=(43520,43520,43520), rgb(DisplayWave1)=(65280,0,0), lsize(DisplayWave1)=2
			ShowInfo
			RenameWindow $S_name, set_cursors
			//Create panel for user definition of UP state 
			NewPanel/k=2/w=(139,341, 422, 522) as "Pause for cursor"
			DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
			AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
			DrawText 60,55,"Scroll through waves"
			DrawText 40,75," to check algorithm performance"
			DrawText 71,134,"then press Continue."
			Button button0,pos={84,140},size={92,20},title="Continue"
			Button button0,proc=UserCursorAdjust_ContButtProc1
			Button button2,pos={90,80},size={30,20},title="<<"
			Button button2,proc=UserScrollBack_ContButtonProc1
			Button button3,pos={150,80},size={30,20},title=">>"
			Button button3,proc=UserScrollFwd_ContButtonProc1
			PauseForUser tmp_PauseforCursor,set_cursors
			killwindow set_cursors
		endif 	
		displayStuff = 0 
		prompt displayStuff, " Happy with result (0) or Run again with different seetings (1)?" 
		doprompt  " inputinput baby" , displayStuff
		if (displayStuff != 1)
			break
			
		endif 
	while(1)
	
	//allow manual check of Up states
	displayStuff = 0 
	prompt displayStuff, " type 1 to manually UPdate UPstate measurements" 
	doprompt  " inputinput baby" , displayStuff
	killwaves/Z filtwav
	if (displayStuff == 1)
		// panda start 
		CurrentWave = 0
		AnalysisWaves1 = wavelist("*DIFUPs*", ";", "") // this your detection signal (0/1) 
		//AnalysisWaves = wavelist("*filt*", ";", "")  // this is your signal
		Wave wav = $StringFromList(CurrentWave, AnalysisWaves)
		Wave UP_marker = $StringFromList(CurrentWave, AnalysisWaves1)
		Duplicate/o $stringfromlist(0,AnalysisWaves) DisplayWave
		duplicate/o $stringfromlist(0,AnalysisWaves1) DisplayWave1	
		Display/w=(220,50, 550, 550) DisplayWave
		ModifyGraph rgb=(0,0,0)
		AppendToGraph/r DisplayWave1
		SetAxis right 0,1
		ShowInfo
		RenameWindow $S_name, set_cursors

		NewPanel/k=2/w=(139,341, 442, 622) as "Pause for cursor"
		DoWindow/C tmp_PauseforCursor        // Set to an unlikely name
		AutoPositionWindow/E/M=1/R=set_cursors    // Put panel near the graph
		DrawText 40,95,"Scroll through waves to check accuracy"
		DrawText 81,174,"then press Finsih."
		Button button0,pos={94,190},size={92,20},title="Finish"
		Button button0,proc=UserCursorAdjust_ContButtProc1
		Button button1,pos={8,108},size={282,20},title="Manually define UP state"
		Button button1,proc=UserEdit_ContButtonProc3
		Button button4,pos={8,128},size={282,20},title="Manually define DOWN state"
		Button button4,proc=UserEdit_ContButtonProc4
		Button button2,pos={90,58},size={30,20},title="<<"
		Button button2,proc=UserScrollBack_ContButtonProc1
		Button button3,pos={170,58},size={30,20},title=">>"
		Button button3,proc=UserScrollFwd_ContButtonProc1
		PauseForUser tmp_PauseforCursor,set_cursors
		DoWindow/k set_cursors
	endif 
	
	//update UP state measurements
	AnalysisWaves1 = wavelist("*DIFUPs*", ";", "") 
	make/o/n=0   UP_sweep 
	Variable baseline1	
	variable counter47 
	for(n = 0; n<ItemsInList(AnalysisWaves1); n+=1)	
		Wave wav = $StringFromList(n, AnalysisWaves)
		Wave UP_marker = $StringFromList(n, AnalysisWaves1)
		make/o/n=(itemsinlist(AnalysisWaves)) UP_frequency
		Make/O/D/N=0 destWave1
		FindLevels/edge=2/q/d=destwave1  UP_marker, 0.5
		make/o/n=10000 UP_amp, UP_duration, UP_area
		duplicate/o wav, filtwav 
		
		Wave W_FindLevels			
		if(V_flag!=2)
			Make/O/D/N=0 coefs;
			FilterIIR/CASC/LO=0.05/ORD=100/COEF coefs, filtwav
			InsertPoints 0,V_LevelsFound, UP_sweep
			UP_sweep[0, (V_LevelsFound-1)] = n
			UP_frequency[n] = V_LevelsFound/60
			counter47  = 0 
			Make/O/D/N=0 destWave
			FindLevels/edge=1/q/d=destWave  UP_marker, 0.5
			do 
				wavestats/q/r=(destWave[counter47], destwave1[counter47]) filtwav
				findlevel/q UP_amp, 0 //find next zero point to insert new value, ugly but whatver
				UP_amp[V_LevelX]  = V_min
				wavestats/q/r=(destWave[counter47]-0.200, destwave[counter47]) filtwav
				UP_amp[V_LevelX]  = V_Avg - UP_amp[V_LevelX]
				UP_duration[V_LevelX] = destwave1[counter47] - destWave[counter47]
				UP_area[V_LevelX]= area( filtwav, destWave[counter47], destWave1[counter47])
				counter47 +=1 
			while (counter47 < V_LevelsFound) 		
		endif
	endfor
	
		
end

Function/D Median(w, x1, x2)	// Returns median value of wave w
	Wave w
	Variable x1, x2	// range of interest

	Variable result

	Duplicate/R=(x1,x2) w, tempMedianWave	// Make a clone of wave
	Sort tempMedianWave, tempMedianWave	// Sort clone
	SetScale/P x 0,1,tempMedianWave
	result = tempMedianWave((numpnts(tempMedianWave)-1)/2)
	KillWaves tempMedianWave
	print result
	return result
End

function MartinFiltFiltIIR (inputWaveString, fs, hiCut,order)
			string inputWaveString
			Variable fs, hiCut, order
			
			duplicate/o $inputWaveString,  $(inputWaveString + "filt")
			wave filtWave = $(inputWaveString+ "filt")
			duplicate/o filtwave, helpWave
			helpWave =filtWave[numpnts(filtWave)-1-p]
			Concatenate/np/o   {filtWave,helpWave}, helpWave2
			duplicate/o helpWave2, helpWave3 
			Concatenate/np/o   {helpWave2,helpWave3}, helpWave 
			Make/O/D/N=0 coefs; DelayUpdate
			FilterIIR/CASC/HI=(hiCut/fs)/ORD=(order)/COEF coefs, helpWave 
			duplicate/o helpWave, helpWave2 
			
			helpWave2 = helpWave[numpnts(helpWave2)-1-p]
			FilterIIR/CASC/HI=(hiCut/fs)/ORD=(order)/COEF coefs, helpWave 
			helpWave = helpWave2[numpnts(helpWave2)-1-p]	
			duplicate/o/r=[(numpnts(filtwave)*2),(numpnts(filtwave)*3)] helpWave, filtWave
			setscale/p x,0,(1/fs), filtwave
			return filtwave
end


Function ExtractUpStateMeasurements (AnalysisWaves,AnalysisWaves1,sF)
	string AnalysisWaves,AnalysisWaves1
	variable sF
	
	// if no analysiswaves have been specified, allow user input 
	if (itemsinlist(AnalysisWaves) ==0)
		string extensionMain = "*conc*"
		string extensionUPs =  "*UPs*"
		prompt extensionMain, "enter unique identifier for original trace" 
		prompt extensionUPs, "enter unique identifier for UPstate detection trace" 
		doprompt "input baby", extensionMain, extensionUPs
		AnalysisWaves = wavelist(extensionMain, ";","")
		AnalysisWaves1 = wavelist(extensionUPs, ";","")
	endif 
	// set up parameters for main routine 
	variable n 
	make/o/n=(itemsinlist(AnalysisWaves)) UP_frequency
	make/o/n=10000 UP_amp, UP_duration, UP_area
	UP_amp = nan
	UP_duration = nan
	UP_area= nan
	Variable baseline1	
	variable counter47 
	variable currMaxlength = 0 
	variable maxval = 0 
	variable totalUPcount = 0 // this is to keep track of total number of UPstates so far 
	duplicate/o  $StringFromList(0, AnalysisWaves) averageUPstate 
	wave avgUP = averageUPstate
	variable pointsToInsert
	avgUP = 0 
	//MAIN ROUTINE 
	// structure: go thorugh sweeps and find UPstates, then go through all the detected UPstates and measure their properties 
	for(n = 0; n<ItemsInList(AnalysisWaves1); n+=1)	
		Wave wav = $StringFromList(n, AnalysisWaves)
		Wave UP_marker = $StringFromList(n, AnalysisWaves1)	
		Make/O/D/N=0 downTrans
		FindLevels/edge=2/q/d=downTrans  UP_marker, 0.5 // destwave1 will contain x values of downstate transitions		
		if(V_flag!=2) //V_flag =  2 wouild mean no level crossings found 
			duplicate/o wav, filtwav 	
			wave filtwav = filtwav
			Make/O/D/N=0 upTrans  
			FindLevels/edge=1/q/d=upTrans  UP_marker, 0.5// destwave will contain X values of  UPstate  transitions
			if (upTrans[0] > downTrans[0]) // if the first downstate transition happens before the first UPstate transition 
				print "WATCH OUT YOUR UPDOWNSTATE TRACE BEGINS WITH A DOWNSTATE TRANSITION, WHICH MAY AFFECT ALGORITHM PERFORMANCE" 
				deletepoints 0,1,downTrans // ignore the first downstate transition 
			endif 
	 		if (downTrans[numpnts(downTrans)-1] < upTrans[numpnts(upTrans)-1])// if the last upstate transition happens after the last downstate transition 
				print "WATCH OUT YOUR UPDOWNSTATE TRACE ENDS WITH AN UPSTATE TRANSITION, WHICH MAY AFFECT ALGORITHM PERFORMANCE" 
				deletepoints numpnts(downTrans)-1,1,downTrans // ignore the last  upstate transition
			endif 
			if (numpnts(downTrans) != numpnts(upTrans)) // this should not happen anymore now 
				print ("The script EXTRACTUPSTATEMEASUREMENTS produces an error for wave: "  +  StringFromList(n, AnalysisWaves))
				return 0 
			endif
			UP_frequency[n] = numpnts(downTrans)/(numpnts(wav)/sF)
			counter47  = 0 	
			
			do // this will go through all the detected UPstates and measure their properties 
				WaveTransform  abs wav
				wave W_Abs = W_Abs
				wavestats/q/r=(upTrans[counter47], downTrans[counter47]) W_Abs
				maxval = V_max
				wavestats/q/r=(upTrans[counter47]-0.3, upTrans[counter47]) wav
				UP_amp[counter47 + totalUPcount]  = maxval - V_Avg
				UP_duration[counter47+ totalUPcount] = downTrans[counter47] - upTrans[counter47]
				UP_area[counter47+ totalUPcount]= area( filtwav, upTrans[counter47], downTrans[counter47])
				duplicate/o/r=(upTrans[counter47] - 1, (numpnts(wav)*deltax(wav))) wav, tmp 
				wave tmp = tmp 
				if (numpnts(tmp) > currMaxlength) 
					 currMaxlength= numpnts(tmp) 	 
				endif
				wavestats/q wav
				tmp-=V_Avg
				setscale/p x,0, deltax(wav), tmp 
				pointsToInsert = numpnts(wav)-numpnts(tmp)
				insertpoints numpnts(tmp), pointsToInsert, tmp 
				avgUP+= tmp
				counter47 +=1 
			while (counter47 < V_LevelsFound) 	
			totalUPcount += counter47 	
		endif
	endfor	
	deletepoints totalUPcount, numpnts(UP_amp)-totalUPcount, UP_amp 
	deletepoints  totalUPcount, numpnts(UP_duration)-totalUPcount, UP_duration
	deletepoints  totalUPcount, numpnts(UP_area)-totalUPcount, UP_area
	deletepoints  currMaxlength, numpnts(avgUP)-currMaxlength, avgUP
	avgUP/=totalUPCount
	wavestats/q UP_duration 
	variable averDur = V_Avg
	display avgUP 
	ModifyGraph rgb=(0,0,0)
	SetDrawEnv arrowfat= 1.5, xcoord= bottom,ycoord= left,linefgc= (65280,0,0),arrow= 1;
	wavestats/q UP_Amp
	variable Avg_Amp = V_avg 
	wavestats/q/r=(0,1) avgUP
	drawline 1 + averDur, V_Avg, 1+ averDur, V_Avg - Avg_Amp
	setaxis bottom, 0, averDur +4
	
	
end