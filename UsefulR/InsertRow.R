insert.row <- function(df, newrow, position) {
  df[seq(position + 1, nrow(df) + 1), ] <- df[seq(position, nrow(df)), ]
  df[position, ] <- newrow
  df
}

out.new.2 <- insert.row(df = out.2, newrow = c("    missing", "0(0%)"), position = 25)
out.new.2 <- insert.row(df = out.new.2, newrow = c("    missing", "0(0%)"), position = 29)
