(ns ribelo.quant
  (:require
   [clojure.test :as t]
   [ribelo.lubeck.quant :as quant]))

(def x [0.003 0.026 0.015 -0.009 0.014 0.024 0.015 0.066 -0.014 0.039])

(t/testing "ribelo.lubeck.quant"
  (t/testing "ann-return-geometric"
    (t/is (= 0.2338146820656939 (quant/ann-return-geometric 12 x))))
  (t/testing "ann-return-simple"
    (t/is (= 0.2148 (quant/ann-return-simple 12 x))))
  (t/testing "annualized-risk"
    (t/is (= 0.08047276972160623 (quant/annualized-risk 12 x))))
  (t/testing "sharpe-ratio"
    (t/is (= 0.7705391416932597  (quant/sharpe-ratio 0 12 x))))
  (t/testing "annualized-sharpe-ratio"
    (t/is (= 2.9055130434129532 (quant/annualized-sharpe-ratio 0 12 x))))
  (t/testing "downside-risk"
    (t/is (= 0.0052630789467763076 (quant/downside-risk 0 x))))
  (t/testing "sortino-ratio"
    (t/is (= 3.4010510161478655 (quant/sortino-ratio 0 0 x))))
  (t/testing "drawdown"
    (t/is (= 0.023000000000000034 (reduce + (quant/drawdown x)))))
  (t/testing "continuous-drawdown"
    (t/is (= 0.02300000000000002 (reduce + (quant/continuous-drawdown x)))))
  (t/testing "average-drawdown"
    (t/is (= 0.01150000000000001 (quant/average-drawdown x))))
  (t/testing "maximum-drawdown"
    (t/is (= 0.014000000000000012 (quant/maximum-drawdown x))))
  (t/testing "rate-of-return"
    (t/is (= 0.18779277315203946 (quant/rate-of-return x))))
  (t/testing "cagr"
    (t/is (= 0.22938756017127182 (quant/cagr (/ 10.0 12.0) x))))
  (t/testing "calmar-ratio"
    (t/is (= 16.70104871897812 (quant/calmar-ratio 0.0 12  x)))))
