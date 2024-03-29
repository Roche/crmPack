
# Stopping the study if the MTD estimation is precise enough, i.e. if robust
# coefficient of variation of the MTD is below 40%.
my_stopping <- StoppingMTDCV(target = 0.3, thresh_cv = 40)
