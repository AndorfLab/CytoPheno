integer_breaks <- function(n) {
  breaker <- scales::pretty_breaks(n)
  function(x) {
    breaks <- breaker(x)
    breaks[breaks == floor(breaks)]
  }
}
