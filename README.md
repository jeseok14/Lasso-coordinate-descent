# Lasso-coordinate-descent
R을 이용한 좌표하강법 기반 Lasso 알고리즘 구현

## 목적
- 좌표하강 알고리즘을 통해 lasso solution을 구함.
- 좌표하강 알고리즘을 R을 통해 직접 구현하고 `glmnet` 패키지의 결과와 비교

## Lasso 모형
- 예측력 향상
- 변수 선택

### 회귀 모형

$y_i = \beta_0 + \sum_{j = 1}^{p}\beta_j x_{ij} + \epsilon_i ( i = 1,2,...,n) \quad \epsilon_i \overset{\mathrm{iid}}{\sim} N(0, \sigma^2)$ 가정

#### 표준화한 변수 사용

- 효과
  - 스케일 차이 감소
  - 모형 안전성 향상
  - 해석의 용이성
  - lasso 벌점에서 표준화된 회귀계수 사용
 
$y_i^s \sim N(0,1)$, $x_{ij}^s \sim N(0,1)$

기존의 $\beta_j$ 대신 $\alpha_j$ 추정

### Lasso 모형

- 효과
  - 예측력 향상
  - 변수선택
 
임의의 $\lambda >0$에 대해서 벌점화 오차제곱합 $Q(\boldsymbol{\alpha})$(목적함수)을 최소화 하는 모수 $\alpha_j$ 구하기

$$Q(\boldsymbol{\alpha}) = \sum_{i = 1}^{n}(y_i^s - \alpha_0 - \alpha_1 x_{i1}^s - \cdots - \alpha_p x_{ip}^s)^2 + \lambda\sum_{j=1}^{p}|\alpha_j|$$

### 좌표하강 알고리즘

목적함수를 최소화하는 모수를 구하기 위해 좌표하강 알고리즘(coordinate descent algorithm: CDA) 사용

- 복수의 모수를 추정하는 경우, 한개씩 개별모수를 순차적으로 업데이트
  - $j$번 째 모수 $\alpha_j$ 추정 시에 $j$를 제외한 모수는 이전 단계 추정치로 고정
- 수렴할 때까지 반복

$$\hat \alpha_j = \arg_{\alpha_j}\min(\alpha_j - t)^2 + \alpha|\alpha_j|$$

여기서 $t = \sum_{i = 1}^{n}(y_i^s - \sum_{\ell \ne j}\hat \alpha_{\ell}x_{i\ell}^s)x_{ij}^s$이다.

이를 만족하는 해는 soft thresholding과 같이 주어진다.

$\hat \alpha_j = soft(t, \lambda/2)$

