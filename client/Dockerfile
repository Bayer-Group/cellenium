
FROM --platform=$TARGETPLATFORM nginx:1.25-alpine3.18 as nginx-version-identifier
RUN  nginx -v 2>&1 | cut -f 2 -d '/' > /etc/nginx/version.txt
RUN cat /etc/nginx/version.txt
RUN echo $TARGETPLATFORM

FROM --platform=$TARGETPLATFORM ubuntu:lunar as builder
ARG TARGETPLATFORM
WORKDIR /app
COPY --from=nginx-version-identifier /etc/nginx/version.txt ./version.txt
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y libpcre3 libpcre3-dev zlib1g zlib1g-dev openssl libssl-dev wget git gcc make libbrotli-dev
RUN export NGINX_VERSION=$(cat version.txt) && wget "https://nginx.org/download/nginx-${NGINX_VERSION}.tar.gz" && tar -zxf nginx-$NGINX_VERSION.tar.gz && git clone https://github.com/google/ngx_brotli && cd ngx_brotli && git submodule update --init && cd ./../nginx-$NGINX_VERSION && ./configure --with-compat --add-dynamic-module=../ngx_brotli && make modules


FROM node:18 as asset-builder
WORKDIR /usr/src/app
COPY package.json yarn.lock ./
RUN yarn
COPY . ./
RUN yarn build

FROM nginx-version-identifier
ARG TARGETPLATFORM
COPY --from=builder /app/nginx-$NGINX_VERSION/objs/ngx_http_brotli_static_module.so /app/nginx-$NGINX_VERSION/objs/ld-linux-aarch64.so* /etc/nginx/modules/
COPY --from=builder /app/nginx-$NGINX_VERSION/objs/ngx_http_brotli_filter_module.so /etc/nginx/modules/
COPY --from=asset-builder /usr/src/app/build /usr/share/nginx/html
RUN if [[ "$TARGETPLATFORM" == "linux/arm64" ]]; then apk update && apk add --no-cache gcompat; fi;
COPY nginx.conf /etc/nginx

CMD ["nginx", "-g", "daemon off;"]


