# Dundun

**Durojaye, C., Fink, L.**, Wald-Fuhrmann, M., Roeske, T. & Larrouy-Maestri, P. (submitted). Perception of Nigerian talking drum performances as speech-like vs. music-like: the role of familiarity and acoustic cues. 

If using anything from this repository, please cite the paper.


___
## Overview

- [Project Description](#Abstract)
- [Classification of dundun performance as speech-like or music-like](#Classification-of-dundun-performance-as-speech-or-music)
- [Predictors of Listener Classification](#Acoustic-features-and-listener-familiarity-predict-classification)
- [Contact](#contact)

___
## Abstract
It seems trivial to identify sound sequences as music or speech, particularly when the sequences come from different sound sources, such as an orchestra and a human voice. Can we also easily distinguish these categories when the sequence comes from the same sound source? On the basis of which acoustic features? We investigated these questions by examining listeners' classification of sound sequences performed by an instrument intertwining both speech and music: the dùndún talking drum. The dùndún is commonly used in south-west Nigeria as a musical instrument but is also perfectly fit for linguistic usage in what has been described as speech surrogates in Africa. One hundred seven participants from diverse geographical locations (15 different mother tongues represented) took part in an online experiment. Fifty-one participants reported being familiar with the dùndún talking drum, 55% of those being speakers of Yorùbá. During the experiment, participants listened to 30 dùndún samples of about 7 seconds long, performed either as music or Yorùbá speech surrogate (n = 15 each) by a professional musician, and were asked to classify each sample as music or speech-like. The classification task revealed the ability of the listeners to identify the samples as intended by the performer, particularly when they were familiar with the dùndún. A logistic regression predicting participants’ classification of the samples from several acoustic features confirmed the perceptual relevance of intensity, timbre, and timing, while showing the interaction of familiarity with features such as pitch. This study provides empirical evidence supporting the discriminating role of acoustic features and the modulatory role of familiarity in teasing apart speech and music.


___
## Classification of dundun performance as speech or music
Participants clearly separated the stimuli into two distinct speech vs. music categories that largely aligned with the intention of the performer. We observed that only four participants categorized every sample as music (solid blue rows near the bottom of the plot), whereas the large majority showed few confusions. Twelve participants exhibited perfect classification (top rows). In A of the figure below, within the speech and music categories, stimuli (columns) are sorted by the number of errors made per stimulus (i.e., the left-most column, stimulus 13M, was least often confused for speech, while 3M was most often confused for speech). Within the speech category, 13S was most clearly perceived as speech, while 5S was most often confused for music. Note that <a href="https://edmond.mpdl.mpg.de/imeji/collection/ovmWl7rLtIiGSv1v" target="_blank">`all dùndún recordings can be accessed online`</a>

A confusion matrix for perceived vs. intended music and speech-like performances is plotted in B below. Overall, the average accuracy of participants on the task was 66%. The average rate people perceived speech when the performance was intended to be music was 12%, while the average rate at which people perceived music when it was intended to be speech was 29%. Collectively, these latter two rates indicate that participants were more likely to perceive speech as music than music as speech. The illustration of confidence ratings (underlying histograms in B, with ratings from 1 to 4) showed similar patterns, with high confidence even in the case of false classifications. Note, however, that the listeners who were unfamiliar (gray) with the dùndún might be less confident when judging speech stimuli (as intended by the performer and confirmed by the three experts).   

Given the imbalance in perceiving speech vs. music, and the statistical properties outlined in the Methods, our main metric of interest for participants’ performance on the task was the Matthews Correlation Coefficient (MCC). An MCC of 1 indicates perfect performance, 0 chance, and -1 perfect misclassification. Participants’ average MCC was 0.61 (+/- 0.33). Participants who were familiar with the dùndún exhibited a significantly higher MCC, compared to those who were unfamiliar with the dùndún (see C below). Also note that participants' familiarity it plotted in the left-hand color bar in A. 

![image](/images/dundun_fig5_rev1.png)

Figure 5. Participants’ classification of music vs. speech-like dùndún performances. (A) Each participants’ (vertical axis) judgment of each stimulus (horizontal axis) as music (blue) or speech (red). Participant performance is sorted descending from least to most errors. Columns are sorted from left to right, within each category (speech and music), from least to most errors per stimulus. The color bar to the left of the plot indicates whether each participant was familiar (golden) or unfamiliar (gray) with the dùndún and whether they spoke Yorùbá (black) or not (white). (B) Confusion matrix for perceived vs. intended stimulus classes, with histograms of participants’ confidence ratings (1-4) for each response type (i.e., quadrant), grouped by familiarity (gold = familiar; gray = unfamiliar). Confidence rating densities were computed within each response type (quadrant) for each familiarity group separately. Means for unfamiliar (gray) and familiar (gold) groups are displayed in the lower left and right corners of each quadrant, respectively. (C) Violin plots and underlying data points indicating the Matthews Correlation Coefficient (MCC) for each participant, separated according to those who were familiar with the dùndún (golden) vs. unfamiliar (gray). The bottom and top horizontal black lines in each distribution represent the 25th and 75th percentiles, the middle line represents the median. An MCC of 1 indicates perfect classification; 0 represents chance, and -1 perfect misclassification. 

___
## Acoustic features and listener familiarity predict classification
In an effort to understand which acoustic features were most relevant in participants’ perception of the dùndún excerpts as music vs. speech-like, we built a linear mixed effects logistic regression model. The binary dependent variable was participants’ perception (speech = 0, music = 1). On the stimulus level, fixed effects included a variety of acoustic measures for each stimulus (related to intensity, pitch, timbre, and timing). Since we observed an effect of familiarity on MCC, with better classification for listeners who were familiar with the dùndún, we included this variable as a fixed effect and in interaction with all acoustic measures. Confidence ratings were also entered as fixed effects. Random intercepts were included for participants and stimuli. 

The figure below shows the odds ratios and confidence intervals for each fixed effect in the model. Fixed effects with an odds ratio less than 1 (red) indicate that a high value on that feature leads to the perception of speech. Odds ratios greater than 1 (blue) indicate that a high value on that feature leads to the perception of music. Overall, the model had a prediction accuracy of 86% and an MCC of 0.71. 

![image](/images/dundun_fig6_rev1.png)

Figure 6. Odds ratios, with confidence interval (CI), for each fixed effect in a logistic regression model predicting participants’ perception of stimuli as music-like (1) or speech-like (0). The vertical line at 1 indicates no effect (i.e., any fixed effect predictor whose odds ratio CI overlaps 1 does not significantly predict participants’ perception). Fixed effects with an odds ratio less than 1 (red) indicate that a high value on that feature leads to the perception of speech. Odds ratios greater than 1 (blue) indicate that a high value on that feature leads to the perception of music. The significance of each fixed effect is indicated with stars (*** p < .001, ** p < .01, * p < .05). Number of observations: 3210. Familiarity was a binary predictor (0 = no; 1 = yes); confidence ranged from 0 to 3; all other variables were continuous. 


___
## Contact
For questions about the analyses in this repository, feel free to <a href="https://lkfink.github.io/" target="_blank">`Dr. Lauren Fink`</a>

