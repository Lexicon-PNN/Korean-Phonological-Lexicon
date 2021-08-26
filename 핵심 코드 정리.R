# <한국어와 영어 기본 단어들의 빈도에 대한 계량적 분석 R 코드>

# 데이터 불러 오기

raw_data = read.csv(file = file.choose())
data = raw_data$data # 빈도 자료를 data라 불리는 데이터 구조에 할당                    
                    
# 정규성 검정

shapiro.test(data) # Shapiro-Wilk 검정

# 두 빈도 분포의 산점도 (예: 빈도 조사 말뭉치 (93)과 세종 말뭉치 (92)) 

x = data_K2 # 빈도 조사 말뭉치 (93)
y = data_K1 # 세종 말뭉치 (93)

plot(x,y, xlab = "빈도 조사 말뭉치 (93)", ylab ="세종 말뭉치 (93)", main="한국어",
     cex.lab=1.8, cex.axis=1.2, cex.main = 2, 
     type = "p", col = "black", pch = 19, lwd = 2)
grid()

# Kendall's tau 검정 (상관관계 분석)

cor.test(data_K2, data_K1, method = "kendall")

# 단어 빈도 빈도순위 정렬 및 빈도 그래프

prop.data = prop.table(data) # 상대빈도
prop.data_rank = rank(-prop.data, ties.method = "first") # 빈도 순위
x = prop.data_rank # x축 빈도 순위
y = prop.data # y축 상대 빈도

plot(x,y, xlab = "빈도 순위", ylab ="상대 빈도", main="BNC (193)", # 빈도 분포 그래프 (BNC (193) 예)
     cex.lab=1.8, cex.axis=1.2, cex.main = 2, 
     type = "p", col = "black", pch = 19, lwd = 2)
grid()

plot(x,y,  log = 'xy', xlab = "빈도 순위", ylab ="빈도", main="BNC (193)", # 로그 스케일 플롯 그래프
     cex.lab=1.8, cex.axis=1.2, cex.main = 2,
     type = "p", col = "black", pch = 19, lwd = 2)

grid()


# 계량적 분석
# 1. poweRlaw 패키지

library(poweRlaw) # poweRlaw 패키지 불러 오기 (poweRlaw를 install 한 후 실행하여야 함)

# 2. 멱함수 분포 Power law distribution
# 자료의 최소 빈도값을 x-min으로 설정할 경우

m_pl_n = displ$new(data)
m_pl_n$setXmin(min(data))
pars <- estimate_pars(m_pl_n)
m_pl_n$setPars(pars) 
bs_p1 = bootstrap_p(m_pl_n, xmins = m_pl_n$xmin, xmax = 1e+06, no_of_sims = 1000, threads = 2) 


# x-min을 추정할 경우
m_pl = displ$new(data)
est = estimate_xmin(m_pl) # xmax를 설정할 경우에는 xmax = 에 자료의 최대 빈도값을 기록함
ls(m_pl)
m_pl$setXmin(est) 
pars <- estimate_pars(m_pl) 
m_pl$setPars(pars) 
m_pl 
bs_p2 = bootstrap_p(m_pl, xmax = 1e+06, no_of_sims = 1000, threads = 2)


# 3. 로그 정규 분포 Lognormal distribution

# 자료의 최소 빈도값을 x-min으로 설정할 경우
m_ln_n = dislnorm$new(data)
m_ln_n$setXmin(min(data))
pars <- estimate_pars(m_ln_n)
m_ln_n$setPars(pars)
m_ln_n 
bs_p1 = bootstrap_p(m_ln_n, xmins = m_ln_n$xmin, xmax = 1e+06, no_of_sims = 1000, threads = 2) 


# x-min을 추정할 경우
m_ln = dislnorm$new(data)
est = estimate_xmin(m_ln) # xmax를 설정할 경우에는 xmax = 에 자료의 최대 빈도값을 기록함
m_ln$setXmin(est) # update the dislnorm object object 
pars <- estimate_pars(m_ln) # pars를 추정할 경우
m_ln$setPars(pars)
m_ln 
bs_p2 = bootstrap_p(m_ln, xmax = 1e+06, no_of_sims = 1000, threads = 2) 


# 4. 지수 분포 Exponential distribution 

# 자료의 최소 빈도값을 x-min으로 설정할 경우
m_exp_n = disexp$new(data)
m_exp_n$setXmin(min(data))
pars <- estimate_pars(m_exp_n)
m_exp_n$setPars(pars)
m_exp_n 
bs_p1 = bootstrap_p(m_exp_n, xmins = m_exp_n$xmin, xmax = 1e+06, no_of_sims = 1000, threads = 2) # Lamda = 0.0005489748, p value = 0.004, gof = 0.3521776
bs_p1

# x-min을 추정할 경우

m_exp = disexp$new(data)
est = estimate_xmin(m_exp) # xmax를 설정할 경우에는 xmax = 에 자료의 최대 빈도값을 기록함
m_exp$setXmin(est) # update the dislnorm object object 
pars <- estimate_pars(m_exp) # pars를 추정할 경우
m_exp$setPars(pars) # pars를 추정할 경우
m_exp 
bs_p2 = bootstrap_p(m_exp, xmax = 1e+06, no_of_sims = 1000, threads = 2)

# 5. 푸아송 분포 Poisson distribution 
# 자료의 최소 빈도값을 x-min으로 설정할 경우

m_pois_n = dispois$new(data)
m_pois_n$setXmin(min(data))
pars <- estimate_pars(m_pois_n)
m_pois_n$setPars(pars)
m_pois_n 
bs_p1 = bootstrap_p(m_pois_n, xmins = m_pois_n$xmin, xmax = 1e+06, no_of_sims = 1000, threads = 2) 

#  x-min을 추정할 경우

m_pois = dispois$new(data)
est = estimate_xmin(m_pois) # xmax를 설정할 경우에는 xmax = 에 자료의 최대 빈도값을 기록함
m_pois$setXmin(est) # update the dislnorm object object 
pars <- estimate_pars(m_pois) # pars를 추정할 경우
m_pois$setPars(pars) # pars를 추정할 경우
m_pois 
bs_p2 = bootstrap_p(m_pois, xmax = 1e+06, no_of_sims = 1000, threads = 2) # Lamda = 99174.17, p value = 0.05, gof = 0.5


# 6. 분포 결과 Plot 그리기 (예 BNC (193))

plot(m_pl_n, xlab = "x (빈도)", ylab ="CDF P(x)", main="BNC (193) 전체",
     cex.lab=1.5, cex.axis=1.2, cex.main = 2,
     pch = 19, lwd = 2)
lines(m_pl_n, col = 2, lwd = 4, lty = 1) # red
lines(m_ln_n, col = 3, lwd = 4, lty = 2) # green
lines(m_exp_n, col = 4, lwd = 4, lty = 4) # blue
lines(m_pois_n, col = 6, lwd = 5, lty = 5) # 청록
legend("bottomleft", lty=c(1, 2, 3, 4), col=c(2, 3, 4, 6), legend=c("멱함수 분포","로그정규 분포", "지수 분포", "푸아송 분포"), 
       title="", bty="n", cex=1.3, lwd = 4)
grid()


plot(m_pl, xlab = "x (빈도)", ylab ="CDF P(x)", main="BNC (193) Xmin 설정",
     cex.lab=1.5, cex.axis=1.2, cex.main = 2,
     pch = 19, lwd = 2)
lines(m_pl, col = 2, lwd = 4, lty = 1) # red
lines(m_ln, col = 3, lwd = 4, lty = 2) # green
lines(m_exp, col = 4, lwd = 4, lty = 4) # blue
lines(m_pois, col = 6, lwd = 5, lty = 5) # 청록
legend("bottomleft", lty=c(1, 2, 3, 4), col=c(2, 3, 4, 6), legend=c("멱함수 분포","로그정규 분포", "지수 분포", "푸아송 분포"), 
       title="", bty="n", cex=1.3, lwd = 4)
grid()


# 우도비 검정 
# (예: 멱함수 분포와 로그 정규 분포를 비교하는데, 멱함수 분포의 X-min을 로그 정규 분포에 사용하여 비교함. 
# 이 경우 R > 0이면 멱함수 분포가 더 적합하다는 결론을 내릴 수 있음)

m_ln$setXmin(m_pl$getXmin())
est = estimate_pars(m_ln)
m_ln$setPars(est)
comp = compare_distributions(m_pl, m_ln) 5712









