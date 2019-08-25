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
        "=": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return x === y });
        },
        "^": function(object1, object2) {
            return map2Scalar(object1, object2, function(x, y) { return x && y ? 1 : 0 });
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

    function isArray(anObject) {
        return Object.prototype.toString.call(anObject) === '[object Array]';
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
        var NUMBERS = /[0-9]+([ \t]+[0-9]+)*/g,
            STRING = /'[^'\n]*'/g,
            BLANK = /[ \t]+/g,
            VARIABLE = /[A-Z]+/g,
            ASSIGN = /←/,
            OUTERPRODUCT = /・./;

        function getVector(vectorString) {
            var splitArray = vectorString.split(/[ \t]/),
                result = map(splitArray, parseFloat);

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

