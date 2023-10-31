FROM node:18 as asset-builder
WORKDIR /usr/src/app
COPY package.json yarn.lock ./
RUN yarn
COPY . ./
RUN yarn build

FROM --platform=$TARGETPLATFORM nginx:1.25-alpine3.18
ARG TARGETPLATFORM
RUN if $TARGETPLATFORM == linux/arm64; then apk update && apk add gcompat; fi;
COPY --from=asset-builder /usr/src/app/build /usr/share/nginx/html
COPY nginx.conf /etc/nginx

CMD ["nginx", "-g", "daemon off;"]


