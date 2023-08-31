#!/usr/bin/env bash
# selection coefficient thresholds
EXP_C1="( S < -0.033 )" 
EXP_C2="(( S > -0.033 ) & ( S < -0.001 ))" 
EXP_C3="( S > -0.001 )" 
# load types
L_MSK="( isHet(GEN[${2}].GT) )"
L_EXPR="((isVariant(GEN[${2}].GT)) & (isHom(GEN[${2}].GT)))"

L1=$(./sh/tally_load_step.sh ${1} "(${EXP_C1} & ${L_MSK} )")
L2=$(./sh/tally_load_step.sh ${1} "(${EXP_C2} & ${L_MSK} )")
L3=$(./sh/tally_load_step.sh ${1} "(${EXP_C3} & ${L_MSK} )")
L4=$(./sh/tally_load_step.sh ${1} "(${EXP_C1} & ${L_EXPR} )")
L5=$(./sh/tally_load_step.sh ${1} "(${EXP_C2} & ${L_EXPR} )")
L6=$(./sh/tally_load_step.sh ${1} "(${EXP_C3} & ${L_EXPR} )")

echo -e "${2}\t${L1}\t${L2}\t${L3}\t${L4}\t${L5}\t${L6}"