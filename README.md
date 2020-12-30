# Dundun

**Durojaye, C., Fink, L.**, Wald-Fuhrmann, M., Roeske, T. & Larrouy-Maestri, P. (submitted). Perception of Nigerian talking drum performances as speech-like vs. music-like: the role of familiarity and acoustic cues. 

If using anything from this repository, please the paper.

Note that code to extract the acoustic features reported in the paper is not provided in this repository. The table of acoustic means required to reproduce all statistical analyses is provided in the Data directory, as are all relevant participant data. 


___
## Overview

- [Project Description](#summary)
- [Listener Classification](#Classification-of-dundun-performance-as-speech-or-music-)
- [Predictors of Listener Classification](#Acoustic-features-and-listener-familiarity-predict-classification)
- [Contact](#contact)

___
## Summary
It seems trivial to identify sound sequences as music or speech, particularly when the sequences come from different sound sources, such as an orchestra vs. a voice. Can we also easily distinguish these categories when the sequence comes from the same sound source? On the basis of which acoustic features? We investigated these questions by examining listeners' classification of sound sequences performed by an instrument intertwining both speech and music: the dùndún talking drum. The dùndún is commonly used in south-west Nigeria as a musical instrument but is also perfectly fit for linguistic usage in what has been described as speech surrogates in Africa. 107 participants from diverse geographical locations, aged from 18 to 75 years (M = 39.22, SD = 15.06), including 51 who were familiar with the dùndún, took part in an online experiment. They listened to 30 dùndún samples of about 7 seconds long performed either as music or speech surrogate (n = 15 each) by a professional musician, and were asked to classify each sample as music or speech-like performance. The classification task revealed the ability of the listeners to identify the samples as intended by the performer, particularly when they were familiar with the dùndún. A logistic regression predicting participants’ classification of the samples from the acoustic features confirmed the perceptual relevance of intensity, timbre, and timing, while showing the interaction of familiarity with features such as pitch. This study provides empirical evidence supporting the discriminating role of acoustic features and the modulatory role of familiarity in teasing apart speech and music.

___
## Classification of dundun performance as speech or music
Participants clearly separated the stimuli into two distinct speech vs. music categories that largely aligned with the intention of the performer. We observed that only four participants categorized every sample as music (solid blue rows near the bottom of the plot), whereas the large majority showed few confusions. Twelve participants exhibited perfect classification (top rows). In A of the figure below, within the speech and music categories, stimuli (columns) are sorted by the number of errors made per stimulus (i.e., the left-most column, stimulus 13M, was least often confused for speech, while 3M was most often confused for speech). Within the speech category, 13S was most clearly perceived as speech, while 5S was most often confused for music. Note that <a href="https://edmond.mpdl.mpg.de/imeji/collection/ovmWl7rLtIiGSv1v" target="_blank">`all dùndún recordings can be accessed online`</a>

A confusion matrix for perceived vs. intended music and speech-like performances is plotted in B below. Overall, the average accuracy of participants on the task was 66%. The average rate people perceived speech when the performance was intended to be music was 12%, while the average rate at which people perceived music when it was intended to be speech was 29%. Collectively, these latter two rates indicate that participants were more likely to perceive speech as music than music as speech. The illustration of confidence ratings (underlying histograms in B, with ratings from 1 to 4) showed similar patterns, with high confidence even in the case of false classifications. Note, however, that the listeners who were unfamiliar (gray) with the dùndún might be less confident when judging speech stimuli (as intended by the performer and confirmed by the three experts).   

Given the imbalance in perceiving speech vs. music, and the statistical properties outlined in the Methods, our main metric of interest for participants’ performance on the task was the Matthews Correlation Coefficient (MCC). An MCC of 1 indicates perfect performance, 0 chance, and -1 perfect misclassification. Participants’ average MCC was 0.61 (+/- 0.33). Participants who were familiar with the dùndún exhibited a significantly higher MCC, compared to those who were unfamiliar with the dùndún (see C below). Also note that participants' familiarity it plotted in the left-hand color bar in A. 

![image](/images/dundun_fig4_perception_violinMCC.png)

___
## Acoustic features and listener familiarity predict classification
In an effort to understand which acoustic features were most relevant in participants’ perception of the dùndún excerpts as music vs. speech-like, we built a linear mixed effects logistic regression model. The binary dependent variable was participants’ perception (speech = 0, music = 1). On the stimulus level, fixed effects included a variety of acoustic measures for each stimulus (related to intensity, pitch, timbre, and timing). Since we observed an effect of familiarity on MCC, with better classification for listeners who were familiar with the dùndún, we included this variable as a fixed effect and in interaction with all acoustic measures. Confidence ratings were also entered as fixed effects. Random intercepts were included for participants and stimuli. 

The figure below shows the odds ratios and confidence intervals for each fixed effect in the model. Fixed effects with an odds ratio less than 1 (red) indicate that a high value on that feature leads to the perception of speech. Odds ratios greater than 1 (blue) indicate that a high value on that feature leads to the perception of music. Overall, the model had a prediction accuracy of 86% and an MCC of 0.71. 

![image](/images/dundun_fig5_perception_logisticMod.png)

Greater pulse clarity predicted perception of music. On the note level, greater mean intensity predicted speech, while greater mean pitch predicted music. In terms of changes between notes, greater changes in pitch predicted speech, while greater changes in timbre predicted music. However, mean pitch and mean pitch change both interacted with participants’ familiarity in the opposite direction of the overall effects, suggesting that participants familiar with the dùndún process its pitch and pitch changes differently.

___
## Contact
For questions about the analyses in this repository, feel free to <a href="https://lkfink.github.io/" target="_blank">`Dr. Lauren Fink`</a>

