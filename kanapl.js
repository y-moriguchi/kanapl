/*
 * KANAPL
 *
 * Copyright (c) 2019 Yuichiro MORIGUCHI
 *
 * This software is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 **/
(function(root) {
    var undef = void 0;

    var monadic = {
        "-": function(array) {
            return map(array, function(x) { return -x; });
        },

        "×": function(array) {
            return map(array, function(x) {
                return x > 0 ? 1 : x < 0 ? -1 : 0;
            });
        },

        "÷": function(array) {
            return map(array, function(x) { return 1 / x; });
        },

        "ι": function(object1) {
            if(isNumber(object1)) {
                return iota(object1, 1, 1);
            } else {
                throw new Error("NOT SUPPORTED");
            }
        },

        "ρ": function(object1) {
            return scalarize(rho(object1));
        },

        "〆": function(object1) {
            return transpose(object1);
        },

        "*": function(array) {
            return map(array, function(x) { return Math.exp(x); });
        },

        "★": function(array) {
            return monadic["*"](array);
        },

        "☆": function(array) {
            return map(array, function(x) {
                return Math.log(x);
            });
        },

        "|": function(array) {
            return map(array, function(x) { return Math.abs(x); });
        },

        "「": function(array) {
            return map(array, function(x) { return Math.floor(x); });
        },

        "」": function(array) {
            return map(array, function(x) { return Math.ceil(x); });
        },

        "〇": function(array) {
            return map(array, function(x) { return Math.PI * x; });
        }
    };

    var dyadic = {
        "+": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return x + y; });
        },

        "-": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return x - y; });
        },

        "×": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return x * y; });
        },

        "÷": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return x / y; });
        },

        "ρ": function(vector1, vector2) {
            function gendim(index, dim) {
                var result = [],
                    resultValue,
                    i,
                    nowIndex = index;

                if(dim >= vector1.length) {
                    return {
                        value: vector2[index],
                        index: index + 1
                    };
                } else {
                    for(i = 0; i < vector1[dim]; i++) {
                        resultValue = gendim(index, dim + 1);
                        index = resultValue.index % vector2.length;
                        result[i] = resultValue.value;
                    }
                    return {
                        value: result,
                        index: index
                    };
                }
            }
            return gendim(0, 0).value;
        },

        "*": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return Math.pow(x, y); });
        },

        "★": function(object1, object2) {
            return dyadic["*"](object1, object2);
        },

        "☆": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) {
                return Math.log(y) / Math.log(x);
            });
        },

        "|": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) {
                if(x === 0) {
                    return y;
                } else if(x > 0 && y >= 0) {
                    return y % x;
                } else if(x > 0 && y <= 0) {
                    return y % x + x;
                } else if(x < 0 && y > 0) {
                    return y % x + x;
                } else {
                    return y % x;
                }
            });
        },

        "「": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) {
                return x > y ? x : y;
            });
        },

        "」": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) {
                return x < y ? x : y;
            });
        },

        "〇": function(anInt, array2) {
            var fn;

            function getFn(anInt) {
                switch(anInt) {
                    case 0:   return function(x) { return Math.sqrt(1 - x * x); };
                    case 1:   return Math.sin;
                    case 2:   return Math.cos;
                    case 3:   return Math.tan;
                    case 4:   return function(x) { return Math.sqrt(1 + x * x); };
                    case 5:   return sinh;
                    case 6:   return cosh;
                    case 7:   return tanh;
                    case -1:  return Math.asin;
                    case -2:  return Math.acos;
                    case -3:  return Math.atan;
                    case -4:  return function(x) { return Math.sqrt(-1 + x * x); };
                    case -5:  return Math.asinh ? Math.asinh : K;
                    case -6:  return Math.acosh ? Math.acosh : K;
                    case -7:  return Math.atanh ? Math.atanh : K;
                    default:  return K;
                }
            }
            fn = getFn(anInt);
            return map(array2, fn);
        },

        "=": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return x === y ? 1 : 0 });
        },

        "≠": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return x !== y ? 1 : 0 });
        },

        "<": function(object1, object2) {
            return mat2Scalar(object1, object2, relCheck(function(x, y) { return x < y }));
        },

        "≦": function(object1, object2) {
            return mat2Scalar(object1, object2, relCheck(function(x, y) { return x <= y }));
        },

        ">": function(object1, object2) {
            return mat2Scalar(object1, object2, relCheck(function(x, y) { return x > y }));
        },

        "≧": function(object1, object2) {
            return mat2Scalar(object1, object2, relCheck(function(x, y) { return x >= y }));
        },

        "^": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return x && y ? 1 : 0 });
        },

        "∧": function(object1, object2) {
            return dyadic["^"];
        },

        "v": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return x || y ? 1 : 0 });
        },

        "∨": function(object1, object2) {
            return dyadic["v"];
        },

        "〆": function(vector1, array1) {
            return substituteArray(vector1, array1);
        }
    };

    var operator = {
        ".": function(fn1, fn2) {
            return function(array1, array2) {
                var result = [],
                    rhoArray1 = rho(array1),
                    rhoArray2 = rho(array2),
                    rhoArray2a = rhoArray2[0];

                function generate1(indicesArray1, rhoArray1, indicesArray2, rhoArray2) {
                    var fold1,
                        i;

                    if(rhoArray1.length > 1) {
                        for(i = 0; i < rhoArray1[0]; i++) {
                            generate1(indicesArray1.concat([i]), rhoArray1.slice(1), indicesArray2, rhoArray2);
                        }
                    } else if(rhoArray2.length > 0) {
                        for(i = 0; i < rhoArray2[0]; i++) {
                            generate1(indicesArray1, rhoArray1, indicesArray2.concat([i]), rhoArray2.slice(1));
                        }
                    } else {
                        fold1 = [];
                        for(i = 0; i < rhoArray2a; i++) {
                            fold1[i] = fn2(getIndex(array1, indicesArray1.concat([i])), getIndex(array2, [i].concat(indicesArray2)));
                        }
                        if(indicesArray1.length + indicesArray2.length > 0) {
                            setIndex(result, indicesArray1.concat(indicesArray2), fold(fold1, fn1));
                        } else {
                            result = fold(fold1, fn1);
                        }
                    }
                }

                generate1([], rhoArray1, [], rhoArray2.slice(1));
                return result;
            };
        }
    };

    function K(x) {
        return x;
    }

    function isNumber(anObject) {
        return typeof anObject === "number";
    }

    function isString(anObject) {
        return typeof anObject === "string";
    }

    function isArray(anObject) {
        return Object.prototype.toString.call(anObject) === '[object Array]';
    }

    function isInteger(anObject) {
        return typeof anObject === 'number' && isFinite(anObject) && Math.floor(anObject) === anObject;
    }

    function sinh(x) {
        var y = Math.exp(x);
        return (y - 1 / y) / 2;
    }

    function cosh(x) {
        var y = Math.exp(x);
        return (y + 1 / y) / 2;
    }

    function tanh(x) {
        var y;
        if(x === Infinity) {
            return 1;
        } else if(x === -Infinity) {
            return -1;
        } else {
            y = Math.exp(x);
            return (y - 1) / (y + 1);
        }
    }

    function deepcopy(anObject) {
        var result;

        function copyAll() {
            var i;

            for(i in anObject) {
                if(anObject.hasOwnProperty(i)) {
                    result[i] = deepcopy(anObject[i]);
                }
            }
        }

        if(isArray(anObject)) {
            result = [];
            copyAll();
            return result;
        } else if(typeof anObject === 'object' && anObject !== null) {
            result = {};
            copyAll();
            return result;
        } else {
            return anObject;
        }
    }

    function relCheck(f) {
        return function(x, y) {
            if(isNumber(x) && isNumber(y)) {
                return f(x, y) ? 1 : 0;
            } else {
                throw new Error("DOMAIN ERROR");
            }
        }
    }

    function scalarize(anArray) {
        if(anArray.length === 1 && !isArray(anArray[0])) {
            return anArray[0];
        } else {
            return anArray;
        }
    }

    function getIndex(array1, indexVector) {
        if(isArray(array1)) {
            return getIndex(array1[indexVector[0]], indexVector.slice(1));
        } else {
            return array1;
        }
    }

    function setIndex(array1, indexVector, value) {
        if(array1[indexVector[0]] === undef) {
            array1[indexVector[0]] = [];
        }
        if(indexVector.length > 1) {
            setIndex(array1[indexVector[0]], indexVector.slice(1), value);
        } else {
            array1[indexVector[0]] = value;
        }
    }

    function transpose(array1) {
        var result = [];

        function walk(array1, indexVector) {
            var i;

            if(isArray(array1)) {
                for(i = 0; i < array1.length; i++) {
                    walk(array1[i], indexVector.concat([i]));
                }
            } else {
                setIndex(result, indexVector.reverse(), array1);
            }
        }

        walk(array1, []);
        return result;
    }

    function substituteArray(vector, array1) {
        var result = [],
            rankVector = rho(array1),
            inIndices = [],
            max;

        function checkVector(vector) {
            var i,
                j,
                max = 0;

            for(i = 0; i < vector.length; i++) {
                if(!isInteger(vector[i])) {
                    throw new Error("DOMAIN ERROR");
                }
                max = max < vector[i] ? vector[i] : max;
            }

            if(max <= 0 || max > vector.length) {
                throw new Error("DOMAIN ERROR");
            }

            outer: for(i = 1; i <= max; i++) {
                for(j = 0; j < vector.length; j++) {
                    if(vector[j] === i) {
                        continue outer;
                    }
                }
                throw new Error("DOMAIN ERROR");
            }
            return max;
        }

        function subst(outIndices) {
            var i,
                j,
                dest;

            if(outIndices.length === max) {
                setIndex(result, outIndices, getIndex(array1, inIndices));
            } else {
                dest = vector.indexOf(outIndices.length + 1);
                for(i = 0; i < rankVector[dest]; i++) {
                    for(j = 0; j < vector.length; j++) {
                        if(vector[j] === outIndices.length + 1) {
                            inIndices[j] = i;
                        }
                    }
                    subst(outIndices.concat([i]));
                }
            }
        }

        if(vector.length !== rankVector.length) {
            throw new Error("LENGTH ERROR");
        }
        max = checkVector(vector);
        subst([]);
        return result;
    }

    function iota(times, start, step) {
        var result = [],
            val = start,
            i;

        start = start === undef ? start : 1;
        step = step === undef ? step : 1;
        for(i = 0; i < times; i++, val += step) {
            result[i] = val;
        }
        return result;
    }

    function rho(anArray) {
        var result = [],
            a1 = anArray,
            i;

        for(i = 0; isArray(a1); i++, a1 = a1[0]) {
            result[i] = a1.length;
        }
        return result;
    }

    function map(anObject, action) {
        var result = [],
            i;

        if(!isArray(anObject)) {
            return action(anObject);
        }

        for(i = 0; i < anObject.length; i++) {
            if(isArray(anObject[i])) {
                result[i] = map(anObject[i], action);
            } else {
                result[i] = action(anObject[i]);
            }
        }
        return result;
    }

    function map2(array1, array2, action) {
        var result = [],
            i;

        if(!isArray(array1) && !isArray(array2)) {
            return action(array1, array2);
        } else if(!(isArray(array1) && isArray(array2))) {
            throw new Error("RANK ERROR");
        } else if(array1.length !== array2.length) {
            throw new Error("LENGTH ERROR");
        } 

        for(i = 0; i < array1.length; i++) {
            if(isArray(array1[i])) {
                result[i] = map2(array1[i], array2[i], action);
            } else {
                result[i] = action(array1[i], array2[i]);
            }
        }
        return result;
    }

    function map2Scalar(object1, object2, action) {
        if(isArray(object1) && isNumber(object2)) {
            return map(object1, function(x) {
                return action(x, object2);
            });
        } else if(isNumber(object1) && isArray(object2)) {
            return map(object2, function(x) {
                return action(object1, x);
            });
        } else {
            return map2(object1, object2, action);
        }
    }

    function outerProduct(fn, array1, array2) {
        function walk(array1) {
            var result = [],
                i;

            if(isArray(array1)) {
                for(i = 0; i < array1.length; i++) {
                    result[i] = walk(array1[i]);
                }
                return result;
            } else {
                return map2Scalar(array1, array2, fn);
            }
        }
        return walk(array1);
    }

    function mapSingle(vector, action) {
        var result = [],
            i;

        for(i = 0; i < vector.length; i++) {
            result[i] = action(vector[i]);
        }
        return result;
    }

    function fold(array1, action) {
        var result = null,
            i;

        for(i = 0; i < array1.length; i++) {
            result = result === null ? array1[i] : action(result, array1[i]);
        }
        return result;
    }

    function parseAPL(program, env) {
        var NUMBERS = /￣?[0-9]+([ \t]+￣?[0-9]+)*/g,
            STRING = /'[^'\n]*'/g,
            BLANK = /[ \t]+/g,
            VARIABLE = /[A-Z]+/g,
            ASSIGN = /←/,
            OUTERPRODUCT = /・./;

        function parseAPLFloat(x) {
            var repl = x;

            repl = repl.replace(/￣/, "-");
            return parseFloat(repl);
        }

        function getVector(vectorString) {
            var splitArray = vectorString.split(/[ \t]/),
                result = map(splitArray, parseAPLFloat);

            if(result.length > 1) {
                return result;
            } else {
                return result[0];
            }
        }

        function getString(aString) {
            var result = [],
                i;

            for(i = 1; i < aString.length - 1; i++) {
                result[i - 1] = aString.charAt(i);
            }
            return result;
        }

        function getVariable(varName) {
            if(!env[varName]) {
                throw new Error("VALUE ERROR");
            }
            return env[varName];
        }

        function skipBlank(index, attr) {
            var result;

            BLANK.lastIndex = index;
            result = BLANK.exec(program);
            if(result && result.index === index) {
                return {
                    lastIndex: BLANK.lastIndex,
                    attr: attr
                };
            } else {
                return {
                    lastIndex: index,
                    attr: attr
                };
            }
        }

        function parseRegex(regex, action, index, attr) {
            var result;

            regex.lastIndex = index;
            result = regex.exec(program);
            if(result && result.index === index) {
                return skipBlank(regex.lastIndex, action(result[0]));
            }
        }

        function walkOperator(index, attr) {
            var op2,
                ch;

            ch = program.charAt(index);
            if(operator[ch]) {
                if(index >= program.length || !(op2 = dyadic[program.charAt(index + 1)])) {
                    throw new Error("SYNTAX ERROR");
                }
                return walkOperator(index + 2, operator[ch](attr, op2));
            } else {
                return {
                    lastIndex: index,
                    attr: attr
                };
            }
        }

        function walkFunction(index, attr) {
            var result,
                resultOp,
                outerOp,
                ch;

            if(index >= program.length) {
                return {
                    lastIndex: index,
                    attr: attr
                };
            }

            ch = program.charAt(index);
            if(ch === ")") {
                return {
                    lastIndex: index + 1,
                    attr: attr
                };
            } else if(ch === "(") {
                return walk(index + 1, attr);
            } else if(dyadic[ch]) {
                resultOp = walkOperator(index + 1, dyadic[ch]);
                result = walk(resultOp.lastIndex, attr);
                return {
                    lastIndex: result.lastIndex,
                    attr: resultOp.attr(attr, result.attr)
                };
            } else if(OUTERPRODUCT.test(program.substring(index, index + 2))) {
                if(!(outerOp = dyadic[program.charAt(index + 2)])) {
                    throw new Error("SYNTAX ERROR");
                }
                result = walk(index + 3, attr);
                return {
                    lastIndex: result.lastIndex,
                    attr: outerProduct(outerOp, attr, result.attr)
                };
            } else {
                throw new Error("SYNTAX ERROR");
            }
        }

        function walkAssign(index, attr) {
            var result,
                rval,
                varName;

            if(!(result = parseRegex(VARIABLE, K, index, attr))) {
                return null;
            }
            varName = result.attr;
            if(!ASSIGN.test(program.charAt(result.lastIndex))) {
                return null;
            }
            rval = walk(result.lastIndex + 1, attr);
            if(rval) {
                env[varName] = rval.attr;
                return skipBlank(rval.lastIndex, rval.attr);
            } else {
                return null;
            }
        }

        function walk(index, attr) {
            var result,
                ch;

            if(index >= program.length) {
                return {
                    lastIndex: index,
                    attr: attr
                };
            }

            ch = program.charAt(index);
            if(ch === "(") {
                result = walk(index + 1, attr);
                return walkFunction(result.lastIndex, result.attr);
            } else if(monadic[ch]) {
                result = walk(index + 1, attr);
                return {
                    lastIndex: result.lastIndex,
                    attr: monadic[ch](result.attr)
                };
            } else if(!!(result = walkAssign(index, attr))) {
                return result;
            } else if(!!(result = parseRegex(NUMBERS, getVector, index, attr))) {
                return walkFunction(result.lastIndex, result.attr);
            } else if(!!(result = parseRegex(STRING, getString, index, attr))) {
                return walkFunction(result.lastIndex, result.attr);
            } else if(!!(result = parseRegex(VARIABLE, getVariable, index, attr))) {
                return walkFunction(result.lastIndex, result.attr);
            } else {
                throw new Error("SYNTAX ERROR");
            }
        }
        return walk(0, "NONE").attr;
    }

    function createEnv(env) {
        var innerEnv = env ? deepcopy(env) : {};

        return function(program) {
            return parseAPL(program, innerEnv);
        };
    }

    if(typeof module !== "undefined" && module.exports) {
        module.exports = createEnv;
    } else {
        root["KANAPL"] = createEnv;
    }
})(this);

