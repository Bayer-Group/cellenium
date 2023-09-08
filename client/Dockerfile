FROM node:16 as asset-builder
WORKDIR /usr/src/app
COPY package.json yarn.lock ./
RUN yarn
COPY . ./
RUN yarn build

FROM --platform=$TARGETPLATFORM ubuntu:lunar as builder
ARG TARGETPLATFORM
WORKDIR /app
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y libpcre3 libpcre3-dev zlib1g zlib1g-dev openssl libssl-dev wget git gcc make libbrotli-dev
RUN wget https://nginx.org/download/nginx-1.25.2.tar.gz && tar -zxf nginx-1.25.2.tar.gz && git clone https://github.com/google/ngx_brotli && cd ngx_brotli && git submodule update --init && cd ./../nginx-1.25.2 && ./configure --with-compat --add-dynamic-module=../ngx_brotli && make modules

FROM --platform=$TARGETPLATFORM nginx:1.25-alpine3.18
ARG TARGETPLATFORM
RUN if $TARGETPLATFORM == linux/arm64; then apk update && apk add gcompat; fi;
COPY --from=builder /app/nginx-1.25.2/objs/ngx_http_brotli_static_module.so /app/nginx-1.25.2/objs/ld-linux-aarch64.so* /etc/nginx/modules/
COPY --from=builder /app/nginx-1.25.2/objs/ngx_http_brotli_filter_module.so /etc/nginx/modules/
COPY --from=asset-builder /usr/src/app/build /usr/share/nginx/html
COPY nginx.conf /etc/nginx

CMD ["nginx", "-g", "daemon off;"]


